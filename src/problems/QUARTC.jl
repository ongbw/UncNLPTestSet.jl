#    Problem : GROUP A
#    *********	 
#    A simple quartic function.
#
#    Origonal SIF Source: problem 157 (p. 87) in
#    A.R. Buckley,
#   "Test functions for unconstrained minimization",
#    TR 1989CS-3, Mathematics, statistics and computing centre,
#    Dalhousie University, Halifax (CDN), 1989.
#
#    QUARTC.SIF classification OUR2-AN-V-0
#
#    Number of variables is variable

function QUARTC_f(x)
	fx = 0.0
	for i in firstindex(x):lastindex(x)
		fx += (x[i] - i)^4
	end
    return fx
end

function QUARTC_g!(x, g)
	for i in firstindex(x):lastindex(x)
		g[i] = 4(x[i] - i)^3
	end
    return g
end

function QUARTC_fg!(x, g)
	fx = 0.0
	for i in firstindex(x):lastindex(x)
		fx  += (x[i] - i)^4
		g[i] = 4(x[i] - i)^3
	end
    return fx, g
end

TestSet["QUARTC"] = UncProgram("QUARTC", QUARTC_f, QUARTC_g!, QUARTC_fg!, 10000, 2ones(10000))