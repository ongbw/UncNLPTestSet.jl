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
∇^2f(\boldmath{x})
```
Determines the full Hessian matrix of an objective function at the point x,
by means of an Automatic Differentiation of the nlp's gradient formula. 
"""
function hessAD(nlp::UncProgram, x::Vector{<:Real}) 
    @lencheck nlp.n x
    g = zeros(length(x))
    return ForwardDiff.jacobian(nlp.g!, g, x)
end


"""
function gAD()

    # TODO: Consider passing Y
"""
function gAD!(nlp, x::Vector{<:Real}, S::Matrix{<:Real}, g::Union{Vector{<:Real}, Nothing}=nothing)
	S_dual = ForwardDiff.Dual{:tag}.(x,  eachcol(S)...)
	# Y_dual = isa(g, Nothing) ? similar(S_dual) : ForwardDiff.Dual{1}.(zeros(nlp.n),  eachcol(S)...) 
	Y_dual = ForwardDiff.Dual{:tag}.(zeros(nlp.n),  eachcol(S)...) 
	# ... when removing g+= statements, we can have Y_dual = similar(S_dual)

	nlp.g!(Y_dual, S_dual)

	Yi = similar(S)
	# update forwardDiff to make extraction more economical
    @views for i in 1:nlp.n
        Yi[i, :] .= Y_dual[i].partials[:]
    end

	if !isa(g, Nothing)
		# Also this should be more economical
		for i in 1:nlp.n
			g[i] = Y_dual[i].value
		end
	end
	return Yi
end

"""
function gAD()

    # bDim corresponds to the number of processors available
    # TODO: we are making the assumption that mod(bDim, nlp.n) ≠ 0, BAD
    # TODO: this should be inplace... 
"""
function gHS(nlp, x, S, bDim::Int)
    nlp.n < bDim && @warn("Block size $bDim ≥ $(nlp.n), the problems dimension")
    bDim = Int(min(nlp.n/2, bDim)) # ensures 2 iterations of gAD to get J(∇f(x))

    # determine the first set of m-1 directions
    m = Int(mod(nlp.n, bDim))
    
    # we determine g on the first iteration
    g = similar(x)
    Y = gAD!(nlp, x, S[:, 1:(m-1)], g)

    # overwrite the last column of S to contain g
    S[:, nlp.n] = g

    for i in m:bDim:(nlp.n-bDim) 
        Yi = gAD!(nlp, x, S[:, i:(i+bDim)])
        Y = [Y Yi] #best way to do this? 
    end
    return g, Y
end

"""
    adjdim!

```math
f(x), ∇f(x)
```
Change the dimensions of a specified problem in the TestSet
"""
function adjdim!(nlp::UncProgram, n::Number=0)
    @warn "This operation may be unstable"
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




export obj, grad, objgrad, adjdim!, hessAD, Programs, SelectProgram, gHS

end # module UncNLPTestSet