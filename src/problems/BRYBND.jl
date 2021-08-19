#    Problem : GROUP A
#    *********
#    The BRYBND problem by Oren.
#
#    Source: problem 31 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
#
#    See also Buckley#73 (p. 41) and Toint#18
#
#    BRYBND.SIF classification SUR2-AN-V-0
#
#    Number of variables is variable

function BRYBND_f(x)
	@warn "Does not agree with CUTEst.jl" 
	mL = 5
	mU = 1
    fx = 0.0
	n = length(x)
    for i in firstindex(x):lastindex(x)-1
		bnd = 0.0
		0 > i-mL ? lo=0 : lo=i-mL 
		n-1 < i+mU ? u=n-1 : u=i+mU

		for j in lo+1:i
        	bnd += x[j]*(1+x[j])
		end
        for j in i+2:u+1
			bnd += x[j]*(1+x[j])
		end

		group = (2.0 + 5.0x[i]^2)x[i] + 1 - bnd
		fx += group^2
    end
    return fx
end

function BRYBND_g!(x, g)
	@warn "Does not agree with CUTEst.jl" 
	mL = 5
	mU = 1
	n = length(x)
    for i in firstindex(x):lastindex(x)-1
		bnd = 0.0
		0 > i-mL ? lo=0 : lo=i-mL 
		n-1 < i+mU ? u=n-1 : u=i+mU

		for j in lo+1:i
        	bnd += x[j]*(1+x[j])
		end
        for j in i+2:u+1
			bnd += x[j]*(1+x[j])
		end

		group = (2.0 + 5.0x[i]^2)x[i] + 1 - bnd
		g[i] += 2.0(2.0 +15.0x[i]^2)group

		for j in lo+1:i
            g[j] -= 2.0(1+ 2.0x[j])group
		end

        for j = i+2:u+1
            g[j] -= 2.0(1+ 2.0*x[j])group
		end
    end
    return g
end

function BRYBND_fg!(x, g)
    @warn "Does not agree with CUTEst.jl" 
	mL = 5
	mU = 1
	n = length(x)
    for i in firstindex(x):lastindex(x)-1
		bnd = 0.0
		0 > i-mL ? lo=0 : lo=i-mL 
		n-1 < i+mU ? u=n-1 : u=i+mU

		for j in lo+1:i
        	bnd += x[j]*(1+x[j])
		end
        for j in i+2:u+1
			bnd += x[j]*(1+x[j])
		end

		group = (2.0 + 5.0x[i]^2)x[i] + 1 - bnd
		fx += group^2
		g[i] += 2.0(2.0 +15.0x[i]^2)group

		for j in lo+1:i
            g[j] -= 2.0(1+ 2.0x[j])group
		end

        for j = i+2:u+1
            g[j] -= 2.0(1+ 2.0*x[j])group
		end
    end
    return fx, g
end

# TestSet["BRYBND"] = UncProgram("BRYBND", BRYBND_f, BRYBND_g!, BRYBND_fg!, 5000, ones(5000))