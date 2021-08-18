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
    item = fx = 0.0
    for i in firstindex(x):5:lastindex(x)
        item = x[i]*i
        fx += item*item

        item = (i+1)*x[i+1]
        fx += item*item

        item = (i+2)*x[i+2]
        fx += item*item

        item = (i+3)*x[i+3]
        fx += item*item

        item = (i+4)*x[i+4]
        fx += item*item
    end
    return fx
end

function POWER_g!(x)
    item = 0.0
    for i in firstindex(x):5:lastindex(x)
        item = x[i]*i
        x[i] = 2.0*item*i

        item = (i+1)*x[i+1]
        x[i+1] = 2.0*item*(i+1)

        item = (i+2)*x[i+2]
        x[i+2] = 2.0*item*(i+2)

        item = (i+3)*x[i+3]
        x[i+3] = 2.0*item*(i+3)

        item = (i+4)*x[i+4]
        x[i+4] = 2.0*item*(i+4)
    end
end

function POWER_fg!(x) 
    item = fx = 0.0
    for i in firstindex(x):5:lastindex(x)
        item = x[i]*i
        fx += item*item
        x[i] = 2.0*item*i

        item = (i+1)*x[i+1]
        fx += item*item
        x[i+1] = 2.0*item*(i+1)

        item = (i+2)*x[i+2]
        fx += item*item
        x[i+2] = 2.0*item*(i+2)

        item = (i+3)*x[i+3]
        fx += item*item
        x[i+3] = 2.0*item*(i+3)

        item = (i+4)*x[i+4]
        fx += item*item
        x[i+4] = 2.0*item*(i+4)
    end
    return fx, x
end

# TODO: Determine an n and a good starting point/min
function POWER(n::Int=25)
    @warn "x0 and minimum not confirmed"
    return UncProgram("POWER", POWER_f, POWER_g!, POWER_fg!, n, zeros(n), zeros(n))
end

export POWER
