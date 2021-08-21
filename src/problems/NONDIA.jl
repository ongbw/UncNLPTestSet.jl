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
#    NONDIA.jl.SIF classification OUR2-AN-V-0
#
#    Number of variables is variable
#
# Daniel Henderson, 08/2021  

function NONDIA_f(x)
	fx = (x[1]-1)^2
	for i in 2:lastindex(x)
		fx += 100*(x[1]-x[i-1]^2)^2
	end
    return fx
end

function NONDIA_g!(x, g)
	g[1] = 2*(x[1]-1)
	for i in 2:lastindex(x)
		g[1]   += 200(x[1]-x[i-1]^2)
		g[i-1] += -400x[i-1]*(x[1] - x[i-1]^2)
	end
    return g
end

function NONDIA_fg!(x, g)
	fx = (x[1]-1)^2
	g[1] = 2*(x[1]-1)
	for i in 2:lastindex(x)
		fx += 100*(x[1]-x[i-1]^2)^2
		g[1]   += 200(x[1]-x[i-1]^2)
		g[i-1] += -400x[i-1]*(x[1] - x[i-1]^2)
	end
    return g
end

TestSet["NONDIA"] = UncProgram("NONDIA", NONDIA_f, NONDIA_g!, NONDIA_fg!, 5000, -ones(5000))