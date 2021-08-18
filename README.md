
# UncNLPTestSet.jl
A testing set for unconstrained non-linear programming problems.  
The testing set consists of a subset of thelarge-unconstrained problems found in the CUTE, _the Constrained and Unconstrained Testing Enviroment_.   

### Under Development
The testing set is currently under development and a user will be exposed to unstable functionality in the main UncNLPTestSet module. However, care has been tacken updating this repository to ensure that unstable features throw appropriate errors, or warnings.  

Once a test program passes adequite testing, it becomes part of the `TestSet`, an exposed dictonary mapping the CUTE problem name to a `UncProgram` instance.
Each problem is named in accordance to their respective .SIF file, where the CUTE problems where origonally encoded.

### Installation
```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/danphenderson/UncNLPTestSet.jl")
julia> using UncNLPTestSet.jl
```

Until a stable version has been released, it is best to regularly update your version of the package after you declaration:

```julia
julia> Pkg.update("UncNLPTestSet.jl")
```

### Interface
During development, it is best to interface with the exported dictionary `TestSet`.
Consider the example of iteratting over the set and calling the core interface at the initial iterate, as specified in the CUTE project. 

```julia
for nlp in Values(TestSet)
    println("Problem name:", nlp.name)
    println("Number of variables:", nlp.n)

    # evaluate the objective function at x0
    obj(nlp, nlp.x0)

    # evaluate the gradient function at x0
    grad(nlp, nlp.x0)

    # evaluate the objective and gradient at x0 over one iteration of the programs dimensions
    obj_grad(nlp, nlp.x0)
end
```

