using UncNLPTestSet, Test, CUTEst, LinearAlgebra, Printf
import NLPModels as jo

for nlp ∈ values(TestSet)
    @printf "\nTesting %s to CUTEst.jl Model at default iterate:\n" nlp.name
    qt_nlp = CUTEstModel(nlp.name) 
    qt_x0  = qt_nlp.meta.x0

    try
        @test norm(obj(nlp, nlp.x0) - jo.obj(qt_nlp, nlp.x0)) < 10e-4
        @test norm(grad(nlp, nlp.x0) - jo.grad(qt_nlp, nlp.x0)) < 10e-4
    catch e
        println(e)
    end

    finalize(qt_nlp)
end
