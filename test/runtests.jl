using UncNLPTestSet, Test, CUTEst, LinearAlgebra
import NLPModels as jo

for nlp ∈ values(TestSet)
    println("\n\nTesting ", nlp.name, " to CUTEst.jl equivalent at their specified x0:")
    qt_nlp = CUTEstModel(nlp.name) 
    qt_x0  = qt_nlp.meta.x0

    try
        @test obj(nlp, nlp.x0) ≈ jo.obj(qt_nlp, qt_x0)
        @test grad(nlp, nlp.x0) ≈ jo.grad(qt_nlp, qt_x0)
    catch e
        println(e)
    end

    finalize(qt_nlp)
end
