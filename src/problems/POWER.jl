#    Problem : GROUP C
#    *********
#    This is an error of a Power problem by Oren.
#
#    Origonal SIF Source:
#    S.S. Oren,
#    Self-scaling variable metric algorithms,
#    Part II: implementation and experiments"
#    Management Science 20(5):863-874, 1974.
#
#    The CUTE implementation specification of the source 
#    formula is not correct. 
#
#    POWER.SIF classification OUR2-AN-V-0
#
# Daniel Henderson, 08/2021  

function POWER_f(x) 
    fx = 0.0
    for i in firstindex(x):lastindex(x)
        fx += x[i]^2*i
    end
    return fx^2
end

function POWER_g!(x, g)
    α = 0.0
    for i in firstindex(x):lastindex(x)
        α += x[i]^2*i
    end
    for i in firstindex(x):lastindex(x)
        g[i] = 4α*i*x[i]
    end
    return g
end

function POWER_fg!(x, g)
    α = 0.0
    for i in firstindex(x):lastindex(x)
        α += x[i]^2*i
    end
    for i in firstindex(x):lastindex(x)
        g[i] = 4α*i*x[i]
    end
    return α^2, g
end

TestSet["POWER"] = UncProgram("POWER", POWER_f, POWER_g!, POWER_fg!, 10000, ones(10000))