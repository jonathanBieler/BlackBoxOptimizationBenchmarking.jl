using BlackBoxOptimizationBenchmarking
using Base.Test

const BBOB = BlackBoxOptimizationBenchmarking

include("../scripts/optimizers_interface.jl")

b = BBOB.benchmark(
    NelderMead(), BlackBoxOptimizationBenchmarking.F1, [100,500,1000], 20, 2, 1e-6,
)

b = BBOB.benchmark(
    NelderMead(), 1, [100,500,1000], 20, 2, 1e-6,
)


@everywhere begin
    using BlackBoxOptimizationBenchmarking
    include(joinpath(Pkg.dir(),"BlackBoxOptimizationBenchmarking/scripts/optimizers_interface.jl"))
end

b = BBOB.benchmark(
    [NelderMead(), GradientDescent()], [1, 2], [100,500,1000], 20, 2, 1e-6,
)

