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
let 

function f(x) 
    fx = 0.0 
    for i in firstindex(x):lastindex(x)-1
        γ = x[i]^2 + x[i+1]^2
        fx += γ^2 + (3.0 - 4.0x[i])
    end
    return fx
end

function g!(x, g)
    for i in firstindex(x):lastindex(x)-1
        γ = x[i]^2 + x[i+1]^2
        g[i] += 4.0γ*x[i] - 4.0
        g[i+1] += 4.0γ*x[i+1]
    end
    return g
end

function fg!(x, g)
    fx = 0.0 
    for i in firstindex(x):lastindex(x)-1
        γ = x[i]^2 + x[i+1]^2
        fx += γ^2 + (3.0 - 4.0x[i])
        g[i] += 4.0γ*x[i] - 4.0
        g[i+1] += 4.0γ*x[i+1]
    end
    return fx, g
end

TestSet["ENGVAL1"] = UncProgram("ENGVAL1", f, g!, fg!, 5000, 2.0*ones(5000))

end # end of local scope