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
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    POWER.SIF classification OUR2-AN-V-0
#
#    Number of variables MUST BE A MULTIPLE OF 5 (due to implemenatation)

function POWER_f(x) 
    fx = 0.0
    for i in firstindex(x):lastindex(x)
        α = x[i]*i
        fx += α^2
    end
    return fx
end

function POWER_g!(x, g)
    for i in firstindex(x):lastindex(x)
        α = x[i]*i
        g[i] = 2.0α
    end
    return g
end

function POWER_fg!(x, g) 
    fx = 0.0
    for i in firstindex(x):lastindex(x)
        α = x[i]*i
        fx += α^2
        g[i] = 2.0α*i
    end
    return fx, g
end

function POWER(n::Int=10000)
    #error("CONFLICT: Gerogian's implemntation does not agree with SIF, but think they are right")
    return UncProgram("POWER", POWER_f, POWER_g!, POWER_fg!, n, ones(n))
end

export POWER
