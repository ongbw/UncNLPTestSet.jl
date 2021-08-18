#    Problem : GROUP A
#    *********
#    This is a simplified version of problem NONDIA.
#
#    Original SIF Source:
#    G. Li,
#    "The secant/finite difference algorithm for solving sparse
#    nonlinear systems of equations",
#    SIAM Journal on Optimization, (to appear), 1990.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    LIARWHD.SIF classification  SUR2-AN-V-0
#
#    Number of variables (at least 2)

function LIARWHD_f(x)  
    fx = 0.0
    for i in firstindex(x):lastindex(x)
        α = 2.0(x[i]^2 - x[1])
        γ = x[i]-1
        fx += α^2 + γ^2
    end
    return fx
end


function LIARWHD_g!(x, g)
    for i in firstindex(x):lastindex(x)
        α = 2.0(x[i]^2 - x[1])
        γ = x[i]-1
        g[i] = 8.0x[i]α + 2.0γ
        g[1] -= 4.0α
    end
    return g
end


function LIARWHD_fg!(x, g)
    fx = 0.0
    for i in firstindex(x):lastindex(x)
        α = 2.0(x[i]^2 - x[1])
        γ = x[i]-1
        fx += α^2 + γ^2
        g[i] = 8.0x[i]α + 2.0γ
        g[1] -= 4.0α
    end
    return fx, g
end

function LIARWHD(n::Int=5000)
    @warn "Minimum not confirmed"
    return UncProgram("LIARWHD", LIARWHD_f, LIARWHD_g!, LIARWHD_fg!, n, 4.0ones(n), ones(n))
end

export LIARWHD