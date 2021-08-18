#   Problem :
#   *********
#
#   The extended Powell singular problem.
#   This problem is a sum of n/4 sets of four terms, each of which is
#   assigned its own group.
#
#   Source:  Problem 13 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
#   See also Toint#19, Buckley#34 (p.85)
#
#   Source: Implementation found here
#   https://www.sfu.ca/~ssurjano/powell.html
#
#   POWELLSG classification: OUR2-AN-V-0
#
#   N is the number of free variables, and should be a multiple of 4

function LANGER_f(x) 
    return error("objective not defined")
end

function LANGER_g!(x)
    return error("gradient not defined")
end

function LANGER_fg!(x) 
    return error("gradient not defined")
end


function LANGER(n::Int=25)
    @warn "x0 and minimum not confirmed"
    return UncProgram("LANGER", LANGER_f, LANGER_g!, LANGER_fg!, n)
end