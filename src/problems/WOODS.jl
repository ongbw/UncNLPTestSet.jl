#    Problem : GROUP A
#    *********
#	 The extended Woods problem.
#	 This problem is a sum of n/4 sets of 6 terms, each of which is
#    assigned its own group.  For a given set i, the groups are
#    A(i), B(i), C(i), D(i), E(i) and F(i). Groups A(i) and C(i) contain 1
#    nonlinear element each, denoted Y(i) and Z(i).
#
#    The problem dimension is defined from the number of these sets.
#    The number of problem variables is then 4 times larger.
#
#	 This version uses a slightly unorthodox expression of Woods
#    function as a sum of squares (see Buckley)
#
#    Origonal SIF Source: problem 14 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
#
#    See also Toint#27, Buckley#17 (p. 101), Conn, Gould, Toint#7
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    WOODS.SIF classification SUR2-AN-V-0
#
#	 This problem is decomposed in n linear groups, the last n-1 of which
#    are 2 x 2 and singular.
#
#    NS is the number of sets (= n/4)

function WOODS_f(x)
	fx = 0.0
	for i in firstindex(x):4:lastindex(x)
		item1  = x[i+1]-x[i]^2
      	item2  = 1-x[i]
      	item3  = x[i+3]-x[i+2]^2
		item4  = 1-x[i+2]
		item5  = x[i+1]+x[i+3]-2.0
		item6  = x[i+1]-x[i+3]
		fx    += 100item1^2 + item2^2 + 90item3^2 + item4^2 + 10item5^2 + 0.1item6^2
	end
    return fx
end

function WOODS_g!(x, g)
	for i in firstindex(x):4:lastindex(x)
		item1  = x[i+1]-x[i]^2
      	item2  = 1-x[i]
      	item3  = x[i+3]-x[i+2]^2
		item4  = 1-x[i+2]
		item5  = x[i+1]+x[i+3]-2.0
		item6  = x[i+1]-x[i+3]
		g[i]   = -400item1*x[i] - 2item2
		g[i+1] = 200item1 + 20item5 + 0.2item6
		g[i+2] = -360item3*x[i+2] - 2item4	
    	g[i+3] = 180item3 + 20item5 - 0.2item6
	end
    return g
end

function WOODS_fg!(x, g)
	fx = 0.0
	for i in firstindex(x):4:lastindex(x)
		item1  = x[i+1]-x[i]^2
      	item2  = 1-x[i]
      	item3  = x[i+3]-x[i+2]^2
		item4  = 1-x[i+2]
		item5  = x[i+1]+x[i+3]-2.0
		item6  = x[i+1]-x[i+3]
		fx    += 100item1^2 + item2^2 + 90item3^2 + item4^2 + 10item5^2 + 0.1item6^2
		g[i]   = -400item1*x[i] - 2item2
		g[i+1] = 200item1 + 20item5 + 0.2item6
		g[i+2] = -360item3*x[i+2] - 2item4	
    	g[i+3] = 180item3 + 20item5 - 0.2item6
	end
    return fx, g
end

@warn "WOODS will break in adjdim!()"
x0 = [j % 2 == 1 ? -3.0 : 1.0 for j in 1:2500]
TestSet["WOODS"] = UncProgram("WOODS", WOODS_f, WOODS_g!, WOODS_fg!, 2500, x0)