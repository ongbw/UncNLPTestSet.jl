### TODO:
#   - Look into problem dimension change (Params*.jl)
#   - use the function dump() to show internals

module UncNLPTestSet
using LinearAlgebra, Printf
using ForwardDiff

# utility
macro lencheck(l, vars...)
    exprs = Expr[]
    for var in vars
      varname = string(var)
      push!(exprs, :(
        if length($(esc(var))) != $(esc(l))
          throw(DimensionError($varname, $(esc(l)), length($(esc(var)))))
        end
      ))
    end
    Expr(:block, exprs...)
end


"""
    UncProgram

A base parent type of each unconstrained non-linear program
"""
mutable struct UncProgram
    name::String
    f::Function
    g!::Function
    fg!::Function
    init::Function

    # init feilds 
    n::Integer
    x0::Vector{T} where T<:Real

    function UncProgram(name, f, g!, fg!, init)
        n, x0 = init()
        new(name, f, g!, fg!, init, n, x0)
    end 
end


"""
    UncNLPTestSet.TestSet

A dictionary mapping the problem name it's corresponding UncProgram.
The dictionary is the users interface to problems contained in UncNLPTestSet.jl
"""
TestSet = Dict{String, UncProgram}()


"""
    obj

```math
f(x)
````
Evaluate the objective function at a point x.
"""
function obj(nlp::UncProgram, x::Vector{<:Real})
    @lencheck nlp.n x
    return nlp.f(x)
end


"""
    grad

```math
∇f(x)
```
Evaluate the gradient of the objective function at a point x.
"""
function grad(nlp::UncProgram, x::Vector{<:Real})
    @lencheck nlp.n x
    g = zeros(length(x))
    return nlp.g!(g, x)
end


"""
    objgrad

```math
f(x), ∇f(x)
```
Evaluate the gradient and it's objective function at a point x, 
by a iterating over the dimesions once.
"""
function objgrad(nlp::UncProgram, x::Vector{<:Real})
    @lencheck nlp.n x
    g = zeros(length(x))
    return nlp.fg!(x, g)
end


"""
    hessAD

```math
∇²f(x)
```
Determines the full Hessian matrix of an objective function at the point x,
by computing the Jacobian the gradient formula. The Jacobian is computed by
a forward automatic differentation from ForwardDiff.jl. 
"""
function hessAD(nlp::UncProgram, x::Vector{<:Real}) 
    @lencheck nlp.n x
    g = zeros(length(x))
    return ForwardDiff.jacobian(nlp.g!, g, x)
end


"""
    adjdim!

```math
f(x), ∇f(x)
```
Change the dimensions of a specified problem in the TestSet
"""
function adjdim!(nlp::UncProgram, n::Number=0)
    @warn "adjdim!: This operation may be unstable"
    n, x0 = nlp.init(n)
    nlp.n = n
    nlp.x0 = x0
end


"""
    Programs

```math
f(x), ∇f(x)
```
A list of unconstrained nonlinear programming problems in the current testing enviroment. 
"""
function Programs()
    for nlp in values(TestSet)
        @printf "%s with dimension %d\n" nlp.name nlp.n
    end
end


"""
    SelectProgram

```math
f(x), ∇f(x)
```
Returns an instance of UncProgram. 
"""
function SelectProgram(key::String)
    if key ∈ keys(TestSet)
        return TestSet[key]
    end
    @warn "The program $s is not in the testing enviroment" 
end

for p in readdir(joinpath(@__DIR__, "problems"))
    include(joinpath("problems", p))
end


""" ----------- Start of solver routines ---------- """


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
    bfgs

```math
∇²f(x)
```
A variant of the Broyden-Fletcher-Goldfard-Shanno inverse Quasi-Newton update extended
to n-directions, such that the multi-secant equations are satisfied.

blah blah blah, list argument requitements here.
"""
function BFGS(H::Matrix{<:Real}, U::Matrix{<:Real}, V::Matrix{<:Real}, ϵ::Float64)
    UTVᵀ = U*pinv(U'V, ϵ)*V'
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








export obj, grad, objgrad, adjdim!, hessAD, Programs, SelectProgram, gHS, BFGS, orth

end # module UncNLPTestSet