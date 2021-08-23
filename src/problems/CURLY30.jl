#    Problem :
#    *********
#    A banded function with semi-bandwidth 30 and
#    negative curvature near the starting point
#
#    Origonal SIF Source: problem 10 in
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

function CURLY30_f(x)
	n = lastindex(x)
    q = similar(x)
	fx = 0.0
	for i in 1:n
        q[i] = sum(x[j] for j = i:min(i + 30, n))
		fx  += q[i]^4 - 20q[i]^2 - 0.1q[i]
	end
    return fx
end

function CURLY30_g!(x, g)
	n = lastindex(x)
    q = similar(x)
	for i in 1:n
        q[i] = sum(x[j] for j = i:min(i + 30, n))
		for j = i:min(i + 30, n)
        	g[j] += 4q[i]^3 - 40q[i] - 0.1
		end
	end
    return g
end

function CURLY30_fg!(x, g)
	fx = 0.0
	n = lastindex(x)
    q = similar(x)
	for i in 1:n
        q[i] = sum(x[j] for j = i:min(i + 30, n))
		fx  += q[i]^4 - 20q[i]^2 - 0.1q[i]
		for j = i:min(i + 30, n)
        	g[j] += 4q[i]^3 - 40q[i] - 0.1
		end
	end
    return fx, g
end

x0 = [0.0001*i/(10000 + 1) for i in 1:10000]
TestSet["CURLY30"] = UncProgram("CURLY30", CURLY30_f, CURLY30_g!, CURLY30_fg!, 10000, x0)