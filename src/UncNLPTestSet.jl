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
    name::AbstractString
    f::Function
    g!::Function
    fg!::Function
    n::Integer
    x0::AbstractVector{<:Real}
end

"""
```math
f(x)
````
Evaluate the objective function at a point x.
"""
function obj(nlp::UncProgram, x::AbstractVector{<:Real})
    return nlp.f(x)
end

"""
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
```math
f(x), ∇f(x)
```
Evaluate the gradient and it's objective function at a point x, by a iterating over the dimesions once.
"""
function obj_grad(nlp::UncProgram, x::AbstractVector{<:Real})
    @lencheck nlp.n x
    g = zeros(length(x))
    return nlp.fg!(x, g)
end

TestSet = Dict{AbstractString, UncProgram}()
# place problems in module
for p in readdir("src/problems")
    include(joinpath("problems", p))
end



export obj, grad, obj_grad, TestSet

end # module UncNLPTestSet