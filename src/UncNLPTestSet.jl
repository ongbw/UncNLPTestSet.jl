"""
	UncNLPTestSet

A C² subset of high-dimensional and nonlinear CUTEr/st problems in native Julia.
Each UncProgram in the UncNLPTestSet specifies an analytical representation of objective f 
and the corresponding gradient ∇f. Additionally, each program has a method for computing the
(f, ∇f) over an economical iteration of the dimensions.

Future development of the UncNLPTestSet should focus on the specification of a standardized
native julia constrained and unconstrained optimization enviroment. UncNLPTestSet.jl contains
a small translation of CUTE, whish is arguably the standard suite for performing numerical 
expeirments in optimization research. A new standard is needed, to test accelerated schemes
utilizing the SIMD parallel nature of Automatic Differentation. The Julia Langague offers
the needed flexibilty in Automatic Differentation computations, primarly through it's
multiple dispatch design which enables operator overloading.

There are many attempts at creating an underlying data structure to hold Standard Form
Programs. The need for a common underlying data structure is to facilitate the definition
of Optimization solving routines. The most supported modeling structure is JuMP, who
acknowledges the diversity and need for a standard program model in the ecosystem.

The UncNLPTestSet remains agnostic to existing Standard Form programming problem
models, as none support flexible automatic differentation. Rather, it serves as the
testing set for a AD based quasi-newton scheme. See qntrhs.jl
"""
module UncNLPTestSet

using LinearAlgebra, Printf, ForwardDiff

"""
    UncProgram

A base parent type of each unconstrained non-linear program
"""
mutable struct UncProgram
    name::String
    f::Function
    g!::Function
    fg!::Function
    init::Function
    n::Integer
    x0::Vector
    function UncProgram(name, f, g!, fg!, init)
        n, x0 = init()
        new(name, f, g!, fg!, init, n, x0)
    end 
end

include("api.jl")
include("qntrhs.jl")

"""
    UncNLPTestSet.TestSet

A dictionary mapping the problem name it's corresponding UncProgram.
"""
TestSet = Dict{String, UncProgram}()
for p in readdir(joinpath(@__DIR__, "programs"))
    include(joinpath("programs", p))
end



export obj, grad, objgrad, adjdim!, hessAD, Programs, SelectProgram, gHS, BFGS, SR1, orth

end