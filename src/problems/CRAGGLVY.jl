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
#    See also Buckley#18
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    CRAGGLVY.SIF classification  OUR2-AY-V-0
#
#    Number of variables (variable, based on imput M, the no. of group sets) 

function CRAGGLVY_f(x) 
    fx = 1.0
    for i in firstindex(x):2:lastindex(x)-3
        item1   = exp(x[i]) - x[i+1]
        item2   = x[i+1] - x[i+2]
        element = x[i+2] - x[i+3]
        item3   = tan(element)+element
        item4   = x[i+3] - 1
        fx     += item1^4 + 100item2^6 + item3^4 + x[i]^8 + item4^2
    end
    return fx
end

function CRAGGLVY_g!(x, g)
    for i in firstindex(x):2:lastindex(x)-3
        item1   = exp(x[i]) - x[i+1]
        item2   = x[i+1] - x[i+2]
        element = x[i+2] - x[i+3]
        item3   = tan(element)+element
        item4   = x[i+3] - 1
        g[i]   += 4*item1^3*exp(x[i]) + 8*x[i]^7
        g[i+1] += -4item1^3 + 600item2^5
        g[i+2] += -600item2^5 + 4item3^3 * (cos(element)^(-2) + 1)
        g[i+3] += -4.0item3^3*(cos(element)^(-2) + 1.0) + 2*item4
    end
    return g
end

function CRAGGLVY_fg!(x, g)
    fx = 1.0 
    for i in firstindex(x):2:lastindex(x)-3
        item1   = exp(x[i]) - x[i+1]
        item2   = x[i+1] - x[i+2]
        element = x[i+2] - x[i+3]
        item3   = tan(element)+element
        item4   = x[i+3] - 1
        fx     += item1^4 + 100item2^6 + item3^4 + x[i]^8 + item4^2
        g[i]   += 4*item1^3*exp(x[i]) + 8*x[i]^7
        g[i+1] += -4item1^3 + 600item2^5
        g[i+2] += -600item2^5 + 4item3^3 * (cos(element)^(-2) + 1)
        g[i+3] += -4.0item3^3*(cos(element)^(-2) + 1.0) + 2*item4
    end
    return fx, g
end



xint = 2*ones(5000)
xint[1] = 1.0
# TestSet["CRAGGLVY"] = UncProgram("CRAGGLVY", CRAGGLVY_f, CRAGGLVY_g!, CRAGGLVY_fg!, 5000, xint)