using Test
using BlackBoxOptimizationBenchmarking, Optim
import BlackBoxOptimizationBenchmarking: minimizer, minimum, optimize
const BBOB = BlackBoxOptimizationBenchmarking

import Base.string

#interface for Optim
pinit(D) = 10*rand(D).-5
optimize(opt::Optim.AbstractOptimizer, f, D, run_length) =
    Optim.optimize(f, pinit(D), opt, Optim.Options(f_calls_limit=run_length, g_tol=1e-12))
    
string(opt::Optim.AbstractOptimizer) = string(typeof(opt).name.name)

test_functions = BBOB.list_functions()

b = BBOB.benchmark(
    NelderMead(), BBOB.sphere, [100, 500, 1000], 20, 2, 1e-6,
)

b = BBOB.benchmark(
    NelderMead(), test_functions[1], [100, 500, 1000], 20, 2, 1e-6,
)

b = BBOB.benchmark(
    GradientDescent(), test_functions[1:2], [100, 500, 1000], 20, 2, 1e-6,
)
