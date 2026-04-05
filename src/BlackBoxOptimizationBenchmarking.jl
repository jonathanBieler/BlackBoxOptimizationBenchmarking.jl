module BlackBoxOptimizationBenchmarking

    using Distributions, Optimization, Optim
    using StaticArrays, Random
    using LinearAlgebra, RecipesBase

    import Optim: minimum, minimizer
    import Base: show

    export BBOBFunction, BenchmarkSetup, bbob_suite
    
    include("BBOBFunction.jl")
    include("benchmark.jl")
    include("plot_benchmark.jl")
    include("plot_functions.jl")

end # module