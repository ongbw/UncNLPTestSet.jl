#    Problem : GROUP A
#    *********
#	 A multi-dimensional variant of HUMPS, a two dimensional function
#    with a lot of humps. The density of humps increases with the
#    parameter ZETA, making the problem more difficult.
#	 
#	 The problem is nonconvex.
#
#    Origonal SIF Source:
#    Ph. Toint, private communication, 1997
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    GENHUMP.SIF classification OUR2-AN-V-0
#
#    Number of variables is variable

function GENHUMP_f(x) 
	fx = 0.0
    for i in firstindex(x):lastindex(x)-1
        item1   = sin(2.0*x[i])
		item2   = sin(2.0*x[i+1])
		fx     += item1^2item2^2 + 0.05(x[i]^2+x[i+1]^2)
    end
    return fx
end

function GENHUMP_g!(x, g)
    for i in firstindex(x):lastindex(x)-1
        item1   = sin(2.0*x[i])
		item2   = sin(2.0*x[i+1])
		g[i]   += 4item1*cos(2x[i])item2^2 + 0.1x[i]
		g[i+1] += 4.0*item2*cos(2.0*x[i+1])*item1*item1 + 0.1*x[i+1]
    end
    return g
end

function GENHUMP_fg!(x, g) 
	fx = 0.0
    for i in firstindex(x):lastindex(x)-1
        item1   = sin(2.0*x[i])
		item2   = sin(2.0*x[i+1])
		fx     += item1^2item2^2 + 0.05(x[i]^2+x[i+1]^2)
		g[i]   += 4item1*cos(2x[i])item2^2 + 0.1x[i]
		g[i+1] += 4.0*item2*cos(2.0*x[i+1])*item1*item1 + 0.1*x[i+1]
    end
    return fx, g
end

x0 = 506ones(5000)
x0[1] = -506 
TestSet["GENHUMPS"] = UncProgram("GENHUMPS", GENHUMP_f, GENHUMP_g!, GENHUMP_fg!, 5000, x0)