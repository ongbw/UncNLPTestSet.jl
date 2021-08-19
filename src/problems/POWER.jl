#    Problem : GROUP A
#    *********
#    The Power problem by Oren.
#
#    Origonal SIF Source:
#    S.S. Oren,
#    Self-scaling variable metric algorithms,
#    Part II: implementation and experiments"
#    Management Science 20(5):863-874, 1974.
#
#    POWER.SIF classification OUR2-AN-V-0
#
#    Number of variables MUST BE A MULTIPLE OF 5 (due to implemenatation)

function POWER_f(x) 
    @warn "CONFLICT: POWER does not agree with CUTEst.jl"
    fx = 0.0
    for i in firstindex(x):lastindex(x)
        fx += (x[i]*i)^2
    end
    return fx
end

function POWER_g!(x, g)
    @warn "CONFLICT: POWER does not agree with CUTEst.jl"
    for i in firstindex(x):lastindex(x)
        g[i] = 2.0i
    end
    return g
end

function POWER_fg!(x, g)
    @warn "CONFLICT: POWER does not agree with CUTEst.jl"
    fx = 0.0
    for i in firstindex(x):lastindex(x)
        fx += (x[i]*i)^2
        g[i] = 2.0i
    end
    return fx, g
end

# TestSet["POWER"] = UncProgram("POWER", POWER_f, POWER_g!, POWER_fg!, 10000, ones(10000))