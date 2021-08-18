#    Problem : GROUP A
#    *********
#    A quartic problem whose Hessian is an arrow-head (downwards) with
#    diagonal central part and border-width of 1.
#
#    Origonal SIF Source:
#    A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
#    "Performance of a multifrontal scheme for partially separable
#    optimization",
#    Report 88/4, Dept of Mathematics, FUNDP (Namur, B), 1988.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    ARWHEAD.SIF classification OUR2-AN-V-0
#
#    Number of variables is variable

function ARWHEAD_f(x) 
    fx = 0.0
    xn_sqr = x[lastindex(x)]^2
    for i in firstindex(x):lastindex(x)-1
        fx += (x[i]^2 + xn_sqr)^2 + 3.0 - 4.0x[i]
    end
    return fx
end

function ARWHEAD_g!(x)
    gn = 0.0
    xn = x[lastindex(x)]
    xn_sqr = xn^2
    for i in firstindex(x):lastindex(x)-1
        α = x[i]^2 + xn_sqr
        x[i] = 4.0x[i]*α - 4.0
        gn += 4.0*xn*α
    end
    x[lastindex(x)] = gn
    return x
end

function ARWHEAD_fg!(x) 
    gn = fx = 0.0
    xn = x[lastindex(x)]
    xn_sqr = xn^2
    for i in firstindex(x):lastindex(x)-1
        α = x[i]^2 + xn_sqr
        fx += (x[i]^2 + xn_sqr)^2 + 3.0 - 4.0*x[i]
        x[i] = 4.0*x[i]*α - 4.0
        gn += 4.0*xn*α
    end
    x[lastindex(x)] = gn
    return x, fx
end

function ARWHEAD(n::Int=5000)
    @warn "minimum not confirmed"
    return UncProgram("ARWHEAD", ARWHEAD_f, ARWHEAD_g!, ARWHEAD_fg!, n, ones(n), zeros(n))
end

export ARWHEAD