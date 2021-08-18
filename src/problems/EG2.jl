#    Problem : GROUP A
#    *********
#    TA simple nonlinear problem given as an example in Section 1.2.4 of
#    the LANCELOT Manual.
#    The problem is non convex and has several local minima.
#
#    Original SIF Source:
#    A.R. Conn, N. Gould and Ph.L. Toint,
#    "LANCELOT, A Fortran Package for Large-Scale Nonlinear Optimization
#    (Release A)"
#    Springer Verlag, 1992.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    EG2.SIF classification OUR2-AN-1000-0
#
#    Number of variables (at least 2)

function EG2_f(x) 
    fx = 0.0 
    for i in firstindex(x):lastindex(x)-1
        item = x[0] + x[i]^2 - 1
        fx += sin(item)
    end
    return fx
end


function EG2_g!(x)
    g = zeros(length(x))
    g[length(x)] = cos(x[n]^2)x[n]
    for i in firstindex(x):lastindex(x)-1
        item = x[0] + x[i]^2 - 1
        g[1] += cos(item)
        g[i] += 2.0cos(item)x[i]
    end
    return g
end


function EG2_fg!(x)
    fx = 0.0 
    g = zeros(length(x))
    g[length(x)] = cos(x[n]^2)x[n]
    for i in firstindex(x):lastindex(x)-1
        item = x[0] + x[i]^2 - 1
        fx += sin(item)
        g[1] += cos(item)
        g[i] += 2.0cos(item)x[i]
    end
    return fx, g
end


function EG2(n::Int=1000)
    @warn "x0 and minimum not confirmed"
    return UncProgram("EG2", EG2_f, EG2_g!, EG2_fg!, n, zeros(n), zeros(n))
end

export EG2