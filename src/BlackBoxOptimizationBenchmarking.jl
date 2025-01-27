module BlackBoxOptimizationBenchmarking

    using Distributions, Memoize, Optimization, Optim
    using LinearAlgebra, RecipesBase

    import Optim: minimum, minimizer
    import Base: show

    export BBOBFunction, BenchmarkSetup
    
    include("BBOBFunction.jl")
    include("benchmark.jl")
    include("plot_benchmark.jl")
    include("plot_functions.jl")

end # module
