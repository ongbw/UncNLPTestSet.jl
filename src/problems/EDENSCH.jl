#    Problem : GROUP A
#    *********
#    The extended Dennis and Schnabel problem, as defined by Li.
#
#    Original SIF Source:
#    G. Li,
#    "The secant/finite difference algorithm for solving sparse
#    nonlinear systems of equations",
#    SIAM Journal on Optimization, (to appear), 1990.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    EDENSCH.SIF classification  OUR2-AN-V-0
#
#    Number of variables (at least 2)

function EDENSCH_f(x)
    fx = 0.0
    for i in firstindex(x):lastindex(x)-1
        item1 = x[i] - 2
        item2 = x[i]x[i+1] - 2x[i+1]
        item3 = x[i+1] + 1
        fx += 16 + item1^4 + item2^2 + item3^2
    end
    return fx
end

function EDENSCH_g!(x, g)
    for i in firstindex(x):lastindex(x)-1
        item1 = x[i] - 2
        item2 = x[i]*x[i+1] - 2x[i+1]
        item3 = x[i+1] + 1
        g[i] += 4item1^3 + 2x[i+1]item2
        g[i+1] += 2item2*(x[i] - 2) + 2item3
    end
    return g
end

function EDENSCH_fg!(x, g)
    fx = 0.0
    for i in firstindex(x):lastindex(x)-1
        item1 = x[i] - 2
        item2 = x[i]*x[i+1] - 2x[i+1]
        item3 = x[i+1] + 1
        fx += 16 + item1^4 + item2^2 + item3^2
        g[i] += 4item1^3 + 2x[i+1]item2
        g[i+1] += 2item2*(x[i] - 2) + 2item3
    end
    return fx, g
end

# TestSet["EDENSCH"] = UncProgram("EDENSCH", EDENSCH_f, EDENSCH_g!, EDENSCH_fg!, 2000, 8ones(2000))