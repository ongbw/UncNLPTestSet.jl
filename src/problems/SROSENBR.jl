#    Problem : GROUP A
#    *********
#	 The separable extension of Rosenbrock's function.
#
#    Origonal SIF Source:problem 21 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
#
#    SROSENBR.SIF classification SUR2-AN-V-0
#
#    Number of variables is variable

function SROSENBR_f(x)
	fx = 0.0
	for i in firstindex(x):2:lastindex(x)
		t1     = 1 - x[i]
		t2     = 10(x[i+1] - x[i]^2)
		fx    += t1^2 + t2^2
	end
    return fx
end

function SROSENBR_g!(x, g)
	for i in firstindex(x):2:lastindex(x)
		t1     = 1 - x[i]
		t2     = 10(x[i+1] - x[i]^2)
		g[i+1] = 20t2
		g[i]   = -2(x[i] * g[i+1] + t1)
	end
    return g
end

function SROSENBR_fg!(x, g)
	fx = 0.0
	for i in firstindex(x):2:lastindex(x)
		t1     = 1 - x[i]
		t2     = 10(x[i+1] - x[i]^2)
		g[i+1] = 20t2
		g[i]   = -2(x[i] * g[i+1] + t1)
		fx    += t1^2 + t2^2
	end
    return fx, g
end

@warn "SROSENBR will break in adjdim!()"
x0 = [j % 2 == 0 ? -1.2 : 1.0 for j in 1:10000]
TestSet["SROSENBR"] = UncProgram("SROSENBR", SROSENBR_f, SROSENBR_g!, SROSENBR_fg!, 10000, x0)