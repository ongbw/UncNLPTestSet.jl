#    Problem : GROUP A
#    *********
#    Another function with nontrivial groups and repetitious elements.
#
#    Original SIF Source:
#    N. Gould, private communication.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    COSINE.SIF classification OUR2-AN-V-0
#
#    Number of variables (variable)
function COSINE_f(x) 
    fx = 0.0
    for i in firstindex(x):lastindex(x)-1
        α = -0.5x[i+1]+x[i]^2
        fx += cos(α)
    end
    return fx
end

# TODO: Late, might be an issue with in-place operation conversion...  
function COSINE_g!(x, g)
    for i in firstindex(x):lastindex(x)-1
        α = -0.5x[i+1]+x[i]^2
        g[i] -= 2.0sin(α)x[i]
        g[i+1] += 0.5sin(α)
    end
    return g
end

# TODO: Late, might be an issue with in-place operation conversion...  
function COSINE_fg!(x, g)
    fx = 0.0 
    for i in firstindex(x):lastindex(x)-1
        α = -0.5x[i+1]+x[i]^2
        fx += cos(α)
        g[i] -= 2.0sin(α)*x[i]
        g[i+1] += 0.5sin(α)
    end
    return fx, g
end


# TODO: Set up initial conditions 
function COSINE(n::Int=10000)
    return UncProgram("COSINE", COSINE_f, COSINE_g!, COSINE_fg!, n, ones(n))
end

nlp = COSINE()
TestSet[nlp.name] = nlp

export COSINE