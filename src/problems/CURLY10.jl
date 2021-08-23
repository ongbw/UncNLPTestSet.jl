#    Problem :
#    *********
#    A banded function with semi-bandwidth 10 and
#    negative curvature near the starting point
#
#    Origonal SIF Source: problem 8 in
#    L. Luksan, C. Matonoha and J. Vlcek
#    Modified CUTE problems for sparse unconstrained optimization,
#    Technical Report 1081,
#    Institute of Computer Science,
#    Academy of Science of the Czech Republic
#
#    CURLY10.SIF classification OUR2-AN-V-0
#
#    Number of variables, n, is variable such that n ≥ 2
#    The initial iterate x0 in the source is not used,
#    rather, x0 is defined as it is the CUTEst enviroment
#
# Daniel Henderson, 08/2021

function CURLY10_f(x)
	n = lastindex(x)
    q = similar(x)
	fx = 0.0
	for i in 1:n
        q[i] = sum(x[j] for j = i:min(i + 10, n))
		fx  += q[i]^4 - 20q[i]^2 - 0.1q[i]
	end
    return fx
end

function CURLY10_g!(g, x)
	n = lastindex(x)
    q = similar(x)
	for i in 1:n
        q[i] = sum(x[j] for j = i:min(i + 10, n))
		for j = i:min(i + 10, n)
        	g[j] += 4q[i]^3 - 40q[i] - 0.1
		end
	end
    return g
end

function CURLY10_fg!(g, x)
	fx = 0.0
	n = lastindex(x)
    q = similar(x)
	for i in 1:n
        q[i] = sum(x[j] for j = i:min(i + 10, n))
		fx  += q[i]^4 - 20q[i]^2 - 0.1q[i]
		for j = i:min(i + 10, n)
			g[j] += 4q[i]^3 - 40q[i] - 0.1
		end
	end
    return fx, g
end

@warn "CURLY* will break in adjdim!() and n≥2 must hold"

x0 = [0.0001*i/(10000 + 1) for i in 1:10000]
TestSet["CURLY10"] = UncProgram("CURLY10", CURLY10_f, CURLY10_g!, CURLY10_fg!, 10000, x0)