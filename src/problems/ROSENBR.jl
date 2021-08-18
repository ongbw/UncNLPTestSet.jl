# defining in accordance to:
#   https://juliadiff.org/ForwardDiff.jl/stable/user/limitations/

# TODO: RENAME TO EXTROSNBR

function ROSENBR_f(x) 
    return 100.0*(x[i+1] - x[i]^2)^2 + (1.0-x[i])^2
end

function ROSENBR_g!(x)
    error("Not implemented")
end

function ROSENBR_fg!(x) 
    error("Not implemented")
end


function ROSENBR()
    return UncProgram("ROSENBR", ROSENBR_f, ROSENBR_g!, ROSENBR_fg!, 2, [-1.2, 1.0,], [1.0, 1.0])
end

export ROSENBR
