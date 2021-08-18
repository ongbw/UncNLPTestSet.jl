#    Problem : GROUP A
#    *********
#    The chained Rosenbrock function as given by Fletcher.
#
#    Original SIF Source:
#    The second problem given by
#    R. Fletcher,
#    "An optimal positive definite update for sparse Hessian matrices"
#    Numerical Analysis report NA/145, University of Dundee, 1992.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    FLETCHR.SIF classification OUR2-AN-V-0
#
#    The Number of variables is variable.

function FLETCHCR_f(x) 
    fx = 0.0
    for i in firstindex(x):lastindex(x)-1
        γ = x[i+1] - x[i] + 1 - x[i]^2
        fx += γ^2
    end
    return fx
end

function FLETCHCR_g!(x, g)
    error("FLETCHR ∇f not implemented")
end
 
function FLETCHCR_fg!(x, g)
    error("FLETCHR ∇f not implemented")
end

function FLETCHCR(n::Int=1000)
    @warn "THERE IS SOME WEIRD STUFF, error in SIF or from Georgians"
    return UncProgram("FLETCHCR", FLETCHCR_f, FLETCHCR_g!, FLETCHCR_fg!, n)
end

# nlp = FLETCHCR()
# TestSet[nlp.name] = nlp

export FLETCHCR