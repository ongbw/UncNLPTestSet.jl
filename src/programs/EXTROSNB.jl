#    Problem : GROUP A
#    *********
#    The extended Rosenbrock function (nonseparable version).
#
#    Original SIF Source: problem 10 in
#    Ph.L. Toint,
#    "Test problems for partially separable optimization and results
#    for the routine PSPMIN",
#    Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
#
#    See also Buckley#116.  Note that MGH#21 is the separable version.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    EXTROSNB.SIF classification  SUR2-AN-V-0
#
#    Number of variables (at least 2)
#
# Daniel Henderson, 08/2021

f = (x) -> begin
    α = x[1] - 1.0
    fx = α^2
    for i in firstindex(x)+1:lastindex(x)
        α = 10*(x[i] - x[i-1]^2)
        fx += α^2
    end
    return fx
end


g! = (g, x) -> begin
    α = x[1] - 1.0
    g[1] = 2α
    for i in firstindex(x)+1:lastindex(x)
        α = 10*(x[i] - x[i-1]^2)
        g[i-1] += -40.0α*x[i-1]
        g[i] += 20.0α
    end
    return g
end


fg! = (g, x) -> begin
    α = x[1] - 1.0
    fx = α^2
    g[1] = 2α
    for i in firstindex(x)+1:lastindex(x)
        α = 10*(x[i] - x[i-1]^2)
        fx += α^2
        g[i-1] += -40.0α*x[i-1]
        g[i] += 20.0α
    end
    return fx, g
end

init = (n::Int=1000) -> begin
	n < 2 && @warn("EXTROSNB: number of variables must be ≥ 2")
	n = max(n, 2)
    
    return n, -1.0*ones(n)
end

TestSet["EXTROSNB"] = UncProgram("EXTROSNB",  f, g!, fg!, init)
