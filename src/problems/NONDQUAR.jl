#    Problem : GROUP A
#    *********
#	 A nondiagonal quartic test problem.
#	 
#	 This problem has an arrow-head type Hessian with a tridiagonal
#    central part and a border of width 1.
#    The Hessian is singular at the solution.
#
#    Origonal SIF Source:
#    A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
#    "Performance of a multi-frontal scheme for partially separable
#    optimization"
#    Report 88/4, Dept of Mathematics, FUNDP (Namur, B), 1988.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    NONDQUAR.jl.SIF classification OUR2-AN-V-0
#
#    Number of variables is variable

function NONDQUAR_f(x)
	n      = length(x)
	fx     = (x[1]-x[2])^2
	item   = x[n-1]+x[n]
	fx    += item^2
	for i in firstindex(x):lastindex(x)-2
      	fx += item^4
	end
    return fx
end

function NONDQUAR_g!(x, g)
	n      = length(x)
	item   = x[0]-x[1]
	g[1]   = 2item
	g[2]   = -2item
	item   = x[n-1]+x[n]
	g[n-1] = g[n] = 2item
	for i in firstindex(x):lastindex(x)-2
		item   = x[i]+x[i+1]+x[n]
      	fx     += item^4
      	g[i]   += 4item^3
      	g[i+1] += 4item^3
      	g[n]   += 4item^3
	end
    return g
end

function NONDQUAR_fg!(x, g)
	n      = length(x)
	item   = x[1]-x[2]
	fx     = item^2
	g[1]   = 2item
	g[2]   = -2item
	item   = x[n-1]+x[n]
	fx    += item^2
	g[n-1] = g[n] = 2item
	for i in firstindex(x):lastindex(x)-2
		item   = x[i]+x[i+1]+x[n]
      	fx     += item^4
      	g[i]   += 4item^3
      	g[i+1] += 4item^3
      	g[n]   += 4item^3
	end
    return fx, g
end

@warn "NONDQUAR will break in adjdim!()"
x0 = [j % 2 == 0 ? -1 : 1 for j in 1:10000]
TestSet["NONDQUAR"] = UncProgram("NONDQUAR", NONDQUAR_f, NONDQUAR_g!, NONDQUAR_fg!, 10000, x0)