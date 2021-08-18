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
    g = copy(x)
    return nlp.g!(g)
end

function grad!(nlp::UncProgram, x::AbstractVector{<:Real})
    @lencheck nlp.n x
    return nlp.g!(x)
end

function obj_grad(nlp::UncProgram, x::AbstractVector{<:Real})
    @lencheck nlp.n x
    # todo, if a convient formula exists call it
    # else we just compute: return obj(nlp, x), grad(nl, x)
    g = copy(x)
    return nlp.fg!(x)
end

function obj_grad!(nlp::UncProgram, x::AbstractVector{<:Real})
    @lencheck nlp.n x
    # todo, if a convient formula exists call it
    # else we just compute: return obj(nlp, x), grad(nl, x)
    return nlp.fg!(x)
end

for p in readdir("src/problems")
    include(joinpath("problems", p))
end

export obj, grad, grad!, obj_grad, obj_grad!

end # module UncNLPTestSet