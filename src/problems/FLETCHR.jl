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
#    The Number of variables is N.

function FLETCHER_f(x) 
    fx = 0.0
    for i in lastindex(x):firstindex(x)-1
        item = x[i+1] - x[i] + 1 - x[i]^2
        fx += item^2
    end
    return 100.0fx
end

# TODO: Late, might be an issue with in-place operation conversion...  
function FLETCHER_g!(x)
    g = zeros(length(x)) 
    for i in lastindex(x):firstindex(x)-1
        item = x[i+1] - x[i] + 1 - x[i]^2
        g[i] += 20.0item*(-2.0x[i] - 1.0)
        g[i+1] += 20.0item
    end
    return 100.0fx, g
end

# TODO: Late, an issue with in-place operation conversion...  
function FLETCHER_fg!(x)
    fx = 0.0
    g = zeros(length(x)) 
    for i in lastindex(x):firstindex(x)-1
        item = x[i+1] - x[i] + 1 - x[i]^2
        fx += item^2
        g[i] += 20.0item*(-2.0x[i] - 1.0)
        g[i+1] += 20.0item
    end
    return 100.0fx, g
end

# TODO: Set up initial conditions 
function FLETCHER(n::Int=10000)
    @warn "x0 and minimum not confirmed"
    return UncProgram("FLETCHR", FLETCHER_f, FLETCHER_g!, FLETCHER_fg!, n, zeros(n), zeros(n))
end

export FLETCHER