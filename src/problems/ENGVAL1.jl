#    Problem : GROUP A
#    *********
#    The ENGVAL1 problem.
#    This problem is a sum of 2n-2 groups, n-1 of which contain 2 nonlinear
#    elements
#
#    Original SIF Source: problem 31 in
#    Ph.L. Toint,
#    "Test problems for partially separable optimization and results
#    for the routine PSPMIN",
#    Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
#
#    See also Buckley#172 (p. 52)
#    SIF input: Ph. Toint and N. Gould, Dec 1989.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    ENGVAL1.SIF classification OUR2-AN-V-0
#
#    N is the number of variables

function ENGVAL1_f(x) 
    fx = 0.0 
    for i in firstindex(x):lastindex(x)-1
        item = x[i]^2 + x[i+1]^2
        fx += item^2 + (3 - 4.0x[i])
    end
    return fx
end


function ENGVAL1_g!(x)
    g = zeros(length(x))
    for i in firstindex(x):lastindex(x)-1
        item = x[i]^2 + x[i+1]^2
        g[i] += 4.0item*x[i] - 4.0
        g[i+1] += 4.0item*x[i+1]
    end
    return g
end

# TODO: Late, might be an issue with in-place operation conversion...  
function ENGVAL1_fg!(x)
    fx = 0.0 
    g = zeros(length(x))
    for i in firstindex(x):lastindex(x)-1
        item = x[i]^2 + x[i+1]^2
        fx += item^2 + (3 - 4.0x[i])
        g[i] += 4.0item*x[i] - 4.0
        g[i+1] += 4.0item*x[i+1]
    end
    return fx, g
end

# TODO
function ENGVAL1(n::Int=5000)
    @warn "x0 and minimum not confirmed"
    return UncProgram("ENGVAL1", ENGVAL1_f, ENGVAL1_g!, ENGVAL1_fg!, n, 2.0*ones(n), zeros(n) )
end

export ENGVAL1