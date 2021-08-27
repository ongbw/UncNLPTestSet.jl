"""
	qntrhs.jl

Implements a Quasi-Newton optimization scheme using a generalized
eigenvalue decomposition to solve the Trust-Region subproblem. The updates
are performed with supplemental partial hessian information obtained through
a forward automatic differentation gradient formula.

The sole purpose of the UncNLPTestSet.jl is to test this scheme.
"""

"""
function gAD()

A subroutine of Algorithm 4.1, as implemented in the paper
"""
function gAD(nlp::UncProgram, x::Vector{<:Real}, S::Matrix{<:Real})
	Sdual = ForwardDiff.Dual{1}.(x,  eachcol(S)...)
	Ydual = ForwardDiff.Dual{1}.(zeros(nlp.n),  eachcol(S)...)
	nlp.g!(Ydual, Sdual)
    
    # extract dual 
    Y = similar(S)
    g = similar(x)
    # TODO: use copyto! here 
    @views for i in 1:length(x)
        Y[i, :] .= Ydual[i].partials[:]
        g[i]     = Ydual[i].value
    end

	return g, Y
end


"""
function gHS()

Algorithm 4.2
"""
function gHS(nlp::UncProgram, x::Vector{<:Real}, S::Matrix)
    b = Int((1+size(S,2))/2)
    g, Y₁ = gAD(nlp, x, S[:, 1:b])
    _, Y₂ = gAD(nlp, x, [S[:,(b+1):(2b-1)] g])

    Y = [Y₁ Y₂[:, 1:b-1]]
    h = Y₂[:, size(Y₂, 2)]
    return g, h, Y
end



"""
    BFGS

```math
∇²f(x)
```
A variant of the Broyden-Fletcher-Goldfard-Shanno inverse Quasi-Newton update extended
to n-directions, such that the multi-secant equations are satisfied.

blah blah blah, list argument requitements here.
"""
function BFGS(H::Symmetric{<:Real}, U::Matrix{<:Real}, V::Matrix{<:Real}, ϵ::Float64)
    UTVᵀ = U*pinv(U'*V, ϵ)*V'
	E = I - UTVᵀ
	return UTVᵀ + E*H*E
end


"""
    SR1

```math
∇²f(x)
```
A variant of SR1, a rank-1 Broyed-class Quasi-Newton update extended,
to n-directions and producing an update of rank b.

blah blah blah, list argument requitements here.
"""
function SR1(H::Matrix{<:Real}, U::Matrix{<:Real}, V::Matrix{<:Real}, ϵ::Float64)
    UlessHV = U-H*V
	T = pinv(Symmetric(UlessHV'V), ϵ)
	return H + UlessHV*T*UlessHV'
end



"""
    orth

Orthoganalizes the space spanned by U with respect to the basis of V
"""
function orth(S::Matrix{<:Real})
	return Matrix(qr(S).Q)
end