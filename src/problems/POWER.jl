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
    γ = fx = 0.0
    for i in firstindex(x):5:lastindex(x)
        γ = x[i]*i
        fx += γ*γ

        γ = (i+1)*x[i+1]
        fx += γ*γ

        γ = (i+2)*x[i+2]
        fx += γ*γ

        γ = (i+3)*x[i+3]
        fx += γ*γ

        γ = (i+4)*x[i+4]
        fx += γ*γ
    end
    return fx
end

function POWER_g!(x, g)
    γ = 0.0
    for i in firstindex(x):5:lastindex(x)
        γ = x[i]*i
        g[i] = 2.0*γ*i

        γ = (i+1)*x[i+1]
        g[i+1] = 2.0*γ*(i+1)

        γ = (i+2)*x[i+2]
        g[i+2] = 2.0*γ*(i+2)

        γ = (i+3)*x[i+3]
        g[i+3] = 2.0*γ*(i+3)

        γ = (i+4)*x[i+4]
        g[i+4] = 2.0*γ*(i+4)
    end
    return g
end

function POWER_fg!(x, g) 
    γ = fx = 0.0
    for i in firstindex(x):5:lastindex(x)
        γ = i*x[i]
        fx += γ^2
        g[i] = 2.0*γ*i

        γ = (i+1)*x[i+1]
        fx += γ^2
        g[i+1] = 2.0*γ*(i+1)

        γ = (i+2)*x[i+2]
        fx += γ^2
        g[i+2] = 2.0*γ*(i+2)

        γ = (i+3)*x[i+3]
        fx += γ^2
        g[i+3] = 2.0*γ*(i+3)

        γ = (i+4)*x[i+4]
        fx += γ^2
        g[i+4] = 2.0*γ*(i+4)
    end
    return fx, g
end

function POWER(n::Int=10000)
    error("The Georgian implementation does not agree with CUTEst")
    @warn "Minimum not confirmed"
    if (k = n % 5) != 0
        n += 5 - k
        @warn("The number of variables given in POWER(n::Int) must be a multiple of 5")
        println("The number of variables has been updated to: $n")
    end
    return UncProgram("POWER", POWER_f, POWER_g!, POWER_fg!, n, ones(n), zeros(n))
end

export POWER
