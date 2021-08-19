### TODO:
#   - Look into problem dimension change (Params*.jl)
#   - use the function dump() to show internals

module UncNLPTestSet
using LinearAlgebra

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
struct UncProgram
    name::String
    f::Function
    g!::Function
    fg!::Function
    n::Int
    x0::AbstractVector{<:Real}
end


"""
    TestSet

A dictionary mapping the problem name it's corresponding UncProgram.
The dictionary is the users interface to problems contained in UncNLPTestSet.jl
"""
TestSet = Dict{AbstractString, UncProgram}()


"""
    obj

```math
f(x)
````
Evaluate the objective function at a point x.
"""
function obj(nlp::UncProgram, x::AbstractVector{<:Real})
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
function grad(nlp::UncProgram, x::AbstractVector{<:Real})
    @lencheck nlp.n x
    g = zeros(length(x))
    return nlp.g!(x, g)
end


"""
    obj_grad

```math
f(x), ∇f(x)
```
Evaluate the gradient and it's objective function at a point x, 
by a iterating over the dimesions once.
"""
function objgrad(nlp::UncProgram, x::AbstractVector{<:Real})
    @lencheck nlp.n x
    g = zeros(length(x))
    return nlp.fg!(x, g)
end

"""
    adjdim!

```math
f(x), ∇f(x)
```
Change the dimensions of an UncProgram in the TestSet
"""
function adjdim!(nlp::UncProgram, n::Number=0)
    n = convert(Int, n)
    if n >= nlp.n 
        x0 = cat(nlp.x0, nlp.x0[nlp.n]ones(n-nlp.n), dims=1) 
        TestSet[nlp.name] = UncProgram(nlp.name, nlp.f, nlp.g!, nlp.fg!, n, x0)
    else # shrink x0
        TestSet[nlp.name] = UncProgram(nlp.name, nlp.f, nlp.g!, nlp.fg!, n, nlp.x0[1:n])
    end
end ## SO FAR stable, but there is a better way to do this. 

for p in readdir(joinpath(@__DIR__, "problems"))
    include(joinpath("problems", p))
end

export obj, grad, obj_grad, TestSet, adjdim!

end # module UncNLPTestSet