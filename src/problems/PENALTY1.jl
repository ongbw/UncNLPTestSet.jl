#    Problem : GROUP A
#    *********	 
#	 This problem is a sum of n+1 least-squares groups, the first n of
#    which have only a linear element.
#    It Hessian matrix is dense.
#
#    Origonal SIF Source: Problem 23 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    PENALTY1.jl.SIF classification SUR2-AN-V-0
#
#    Number of variables is variable

function PENALTY1_f(x)
	a  = 1e-5
	fx = tail = 0.0
	for i in firstindex(x):lastindex(x)
		item  = x[i] - 1
      	fx   += item^2*a
      	tail += x[i]^2
	end
	fx += (tail - 0.25)^2
    return fx
end

function PENALTY1_g!(x, g)
	a  = 1e-5
	tail = 0.0
	for i in firstindex(x):lastindex(x)
		item  = x[i] - 1
      	g[i]  = 2a*item
      	tail += x[i]^2
	end
	for i in firstindex(x):lastindex(x)
		g[i] += 4*(tail-0.25)x[i]
	end
    return g
end

function PENALTY1_fg!(x, g)
	a  = 1e-5
	tail = 0.0
	for i in firstindex(x):lastindex(x)
		item  = x[i] - 1
      	fx   += item^2*a
      	g[i]  = 2a*item
      	tail += x[i]^2
	end
	fx += (tail - 0.25)^2
	for i in firstindex(x):lastindex(x)
		g[i] += 4*(tail-0.25)x[i]
	end
    return fx, g
end

@warn "PENALTY1 will break in adjdim!()"
TestSet["PENALTY1"] = UncProgram("PENALTY1", PENALTY1_f, PENALTY1_g!, PENALTY1_fg!, 1000, Vector(1.0:1000.0))