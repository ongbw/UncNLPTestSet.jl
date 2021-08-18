module UncNLPTestSet

using LinearAlgebra

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
    min::AbstractVector{<:Real}
end

function obj(nlp::UncProgram, x::AbstractVector{<:Real})
    return nlp.f(x)
end

function grad(nlp::UncProgram, x::AbstractVector{<:Real})
    @lencheck nlp.n x
    g = zeros(length(x))
    return nlp.g!(x, g)
end

function obj_grad(nlp::UncProgram, x::AbstractVector{<:Real})
    @lencheck nlp.n x
    g = zeros(length(x)) # similar does not work?
    return nlp.fg!(x, g)
end

for p in readdir("src/problems")
    include(joinpath("problems", p))
end

export obj, grad, obj_grad

end # module UncNLPTestSet