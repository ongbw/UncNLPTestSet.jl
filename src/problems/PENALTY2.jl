#    Problem : GROUP A
#    *********	 
#    The second penalty function
#
#    This is a nonlinear least-squares problem with M=2*N groups.
#    Group 1 is linear.
#    Groups 2 to N use 2 nonlinear elements.
#    Groups N+1 to M-1 use 1 nonlinear element.
#    Group M uses N nonlinear elements.
#    The Hessian matrix is dense.

#
#    Origonal SIF Source: Problem 24 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
#
#    See also Buckley#112 (p. 80)
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    PENALTY2.SIF classification SUR2-AN-V-0
#
#    Number of variables is variable

function PENALTY2_f(x)
	a 	  = 1e-5
	fx    = (x[1]-0.2)^2
	tail  = x[1]^2*(n-1)
	for i in firstindex(x)+1:lastindex(x)
		item1 = exp(x[i]/10) + exp(x[i-1]/10) - exp(i/10.0)-exp((i-1)/10.0)
		item2 = exp(x[i]/10) -exp(-1/10.0)
		fx   += (item1^2 + item2^2)a
		tail   += (n-i-1)x[i]^2
	end
	fx += (tail - 0.25)^4
    return fx
end

function PENALTY2_g!(x, g)
	n     = length(x)
	a 	  = 1e-5
	item1 = x[1]-0.2
	g[1]  = 2item1
	tail += x[1]^2*n
	for i in firstindex(x)+1:lastindex(x)
		item1 = exp(x[i]/10) + exp(x[i-1]/10) - exp(i/10.0)-exp((i-1)/10.0)
		item2   = exp(x[i]/10) -exp(-1/10.0)
		tail   += (n-i-1)x[i]^2
		g[i-1] += 0.2a*exp(x[i-1]/10)item1
		g[i]   += 0.2a*exp(x[i]/10)(item1 + item2)
	end
	for i in firstindex(x):lastindex(x)
		g[i] += 4(tail - 0.25)x[i]
	end
    return g
end

function PENALTY2_fg!(x, g)
	n     = length(x)
	a 	  = 1e-5
	item1 = x[1]-0.2
	fx    = item1^2
	g[1]  = 2item1
	tail += x[1]^2*(n-1)
	for i in firstindex(x)+1:lastindex(x)
		item1 = exp(x[i]/10) + exp(x[i-1]/10) - exp(i/10.0)-exp((i-1)/10.0)
		item2   = exp(x[i]/10) -exp(-1/10.0)
		tail   += (n-i-1)x[i]^2
		fx     += (item1^2 + item2^2)a
		g[i-1] += 0.2a*exp(x[i-1]/10)item1
		g[i]   += 0.2a*exp(x[i]/10)(item1 + item2)
	end
	fx += (tail - 0.25)^4
	for i in firstindex(x):lastindex(x)
		g[i] += 4(tail - 0.25)x[i]
	end
    return fx, g
end


TestSet["PENALTY2"] = UncProgram("PENALTY2", PENALTY2_f, PENALTY2_g!, PENALTY2_fg!, 100, 0.5ones(100))