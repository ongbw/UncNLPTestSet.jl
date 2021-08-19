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
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    NONDIA.jl.SIF classification OUR2-AN-V-0
#
#    Number of variables is variable

function NONDIA_f(x)
	item = x[1]-1
	fx = item^2
	item = 10(x[1] - x[1]^2)
	fx += item^2
	for i in firstindex(x)+1:lastindex(x)
		item = 10(x[1] - x[i]^2)
		fx += item^2
	end
    return fx
end

function NONDIA_g!(x, g)
	g[1] = 2x[1]-2
	item = 10(x[1] - x[1]^2)
	g[1] += (20 - 40x[1])item
	for i in firstindex(x)+1:lastindex(x)
		item = 10(x[1] - x[i]^2)
		g[1] += 20item
		g[i] = -40x[i]item
	end
    return g
end

function NONDIA_fg!(x, g)
	item = x[1]-1
	fx = item^2
	g[1] = 2item
	item = 10(x[1] - x[1]^2)
	fx += item^2
	g[1] += (20 - 40x[1])item
	for i in firstindex(x)+1:lastindex(x)
		item = 10(x[1] - x[i]^2)
		fx += item^2
		g[1] += 20item
		g[i] = -40x[i]item
	end
    return fx, g
end


TestSet["NONDIA"] = UncProgram("NONDIA", NONDIA_f, NONDIA_g!, NONDIA_fg!, 5000, -ones(5000))