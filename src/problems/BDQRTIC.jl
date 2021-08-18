#    Problem : GROUP A
#    *********
#    This problem is quartic and has a banded Hessian with bandwidth = 9.
#
#    Origonal SIF Source:
#    Problem 61 in
#    A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
#    "Performance of a multifrontal scheme for partially separable
#    optimization",
#    Report 88/4, Dept of Mathematics, FUNDP (Namur, B), 1988.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    BDQRTIC.SIF classification SUR2-AN-V-0
#
#    Number of variables (variable)

function BDQRTIC_f(x) 
    fx = 0.0
    xn_sqr = x[lastindex(x)]^2
    for i in firstindex(x):lastindex(x)-4
        α = 3.0 - 4.0x[i]
        γ = x[i]^2 + 2.0x[i+1]^2 + 3.0x[i+2]^2 + 4.0x[i+3]^2 + 5.0xn_sqr
        fx += α*γ + γ^2
    end
    return fx
end

# TODO: This is an instance where an in-place operation would be very bad..
#       fix current implementation when it becomes more clear..  
function BDQRTIC_g!(x)
    g = zeros(length(x))
    xn = x[lastindex(x)]
    xn_sqr = xn^2
    gn = 0.0
    for i in firstindex(x):lastindex(x)-4
        α = 3.0 - 4.0x[i]
        γ = x[i]^2 + 2.0x[i+1]^2 + 3.0x[i+2]^2 + 4.0x[i+3]^2 + 5.0xn_sqr
        g[i] += -8.0α + 4.0γ*x[i];
        g[i+1] += 8.0γ*x[i+1];
        g[i+2] += 12.0γ*x[i+2];
        g[i+3] += 16.0γ*x[i+3];
        gn += 20.0γ*xn;
    end
    g[lastindex(x)] = gn
    return g
end

# TODO: This is an instance where an in-place operation would be very bad..
#       fix current implementation when it becomes more clear..  
function BDQRTIC_fg!(x)
    fx = 0.0 
    g = zeros(length(x))
    xn = x[lastindex(x)]
    xn_sqr = xn^2
    gn = 0.0
    for i in firstindex(x):lastindex(x)-4
        α = 3.0 - 4.0x[i]
        γ = x[i]^2 + 2.0x[i+1]^2 + 3.0x[i+2]^2 + 4.0x[i+3]^2 + 5.0xn_sqr
        fx += α*γ + γ^2
        g[i] += -8.0α + 4.0γ*x[i];
        g[i+1] += 8.0γ*x[i+1];
        g[i+2] += 12.0γ*x[i+2];
        g[i+3] += 16.0γ*x[i+3];
        gn += 20.0γ*xn;
    end
    g[lastindex(x)] = gn
    return g, fx
end


function BDQRTIC(n::Int=5000)
    @warn "minimum not confirmed & Instance of a grad call making a redundant copy of x" 
    return UncProgram("BDQRTIC", BDQRTIC_f, BDQRTIC_g!, BDQRTIC_fg!, n, ones(n), zeros(n))
end

export BDQRTIC