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
    @warn "DOES NOT AGREE WITH CUTEst.jl, but this should be correct"  
    fx = 0.0
    for i in firstindex(x):lastindex(x)-1
        γ = x[i+1] - x[i] + 1 - x[i]^2
        fx += 100*γ^2
    end
    return fx
end

function FLETCHCR_g!(x, g)
    @warn "DOES NOT AGREE WITH CUTEst.jl, but this should be correct" 
    for i in firstindex(x):lastindex(x)-1
        g[i] += 200*(1 + 2x[i])(-1+x[i]+x[i]^2-x[i+1])
        g[i+1] += 200*(1 - x[i] - x[i]^2 + x[i+1])
    end
    return g
end
 
function FLETCHCR_fg!(x, g)
    fx = 0.0
    @warn "DOES NOT AGREE WITH CUTEst.jl, but this should be correct" 
    for i in firstindex(x):lastindex(x)-1
        γ = x[i+1] - x[i] + 1 - x[i]^2
        fx += 100*γ^2
        g[i] += 200*(-1 - 2x[i])*(1 - x[i] - x[i]^2 + x[1+i])
        g[i+1] += 200*(1 - x[i] - x[i]^2 + x[1+i])
    end
    return fx, g
end

# TestSet["FLETCHCR"] = UncProgram("FLETCHCR", FLETCHCR_f, FLETCHCR_g!, FLETCHCR_fg!, 1000, zeros(1000))