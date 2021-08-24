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

    # TODO: this is ridiculous implementation
"""
function gAD!(nlp::UncProgram, x::Vector{<:Real}, S::Matrix{<})
    @lencheck nlp.n size(S, 1)
    dual = ForwardDiff.Dual{1}.(x,  eachcol(S)...)
    result = ForwardDiff.Dual{1}.(zeros(nlp.n),  eachcol(S)...) 
	nlp.g!(result, dual)
    g = similar(x)
    Y = similar(S)
    for i in 1:nlp.n
        Y[i, :] = result[i].partials[:]
        g[i] = result[i].values
    end
    return g, S
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
    select

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




export obj, grad, objgrad, adjdim!, hessAD, Programs, SelectProgram

end # module UncNLPTestSet