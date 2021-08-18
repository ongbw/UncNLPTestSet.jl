#    Problem : GROUP A
#    *********
#    Extended Cragg and Levy problem.
#    This problem is a sum of m  sets of 5 groups,
#    There are 2m+2 variables. The Hessian matrix is 7-diagonal.
#
#    Original SIF Source:
#    Ph. L. Toint,
#    "Test problems for partially separable optimization and results for the routine
#    PSPMIN", Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    CRAGGLVY.SIF classification  OUR2-AY-V-0
#
#    Number of variables (variable, based on imput M, the no. of group sets) 

function CRAGGLVY_f(x) 
    fx = 0.0
    for i in firstindex(x):lastindex(x)-1
        item1 = exp(x[i]) - x[i+1]
        item2 = x[i+1] - x[i+2]
        element = x[i+2] - x[i+3]
        item3 = tan(element)+element
        item4 = x[i+3] - 1
        fx += item1^4 + 100.0item2^6 + item3^4 + x[i]^8 + item4^2
    end
    return fx
end

function CRAGGLVY_g!(x)
    return error("gradient not implented")
end

# TODO: Late, might be an issue with in-place operation conversion...  
function CRAGGLVY_fg!(x)
    return error("gradient not implented")
    fx = 0.0 
    g = zeros(length(x))
    for i in firstindex(x):lastindex(x)-1
        item1 = exp(x[i]) - x[i+1]
        item2 = x[i+1] - x[i+2]
        element = x[i+2] - x[i+3]
        item3 = tan(element)+element
        item4 = x[i+3] - 1
        fx += item1^4 + 100.0item2^6 + item3^4 + x[i]^8 + item4^2
        g[i] += 4.0item1^3exp(x[i]) + 8.0x[i]^7
        g[i+1] += -4.0item1^3 + 600.0item2^5
        g[i+2] += -600.0item2^5 + 4.0item3^3 * (cos(element)^(-2) + 1.0)
        # g[i+3] += -4.0item3^3*(cos(element)^(-2) + 1.0)c+ 2*item4 figure out what the 'c' is 
    end
    return fx, g
end



function CRAGGLVY(n::Int=5000)
    error("Issues")
    xint = 2*ones(5000)
    xint[1] = 1.0
    @warn "x0 and minimum not confirmed"
    return UncProgram("CRAGGLVY", CRAGGLVY_f, CRAGGLVY_g!, CRAGGLVY_fg!, n, xint, zeros(n) ) # NOT CONFIRMED
end

export CRAGGLVY