module BlackBoxOptimizationBenchmarking

    using Distributions, Memoize, Optim
    using LinearAlgebra, RecipesBase

    import Base: show
    import Optim: minimizer, minimum, optimize
    export BBOBFunction

    
    include("BBOBFunction.jl")
    include("benchmark.jl")
    include("plot_benchmark.jl")

end # module
