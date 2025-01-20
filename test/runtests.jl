using Test
using BlackBoxOptimizationBenchmarking, Optim, Plots
import BlackBoxOptimizationBenchmarking: minimizer, minimum, optimize
const BBOB = BlackBoxOptimizationBenchmarking

import Base.string

map(BBOB.test_x_opt, BBOB.list_functions())

#interface for Optim
pinit(D) = 10*rand(D).-5
optimize(opt::Optim.AbstractOptimizer, f, D, run_length) =
    Optim.optimize(f, pinit(D), opt, Optim.Options(f_calls_limit=run_length, g_tol=1e-12))
    
string(opt::Optim.AbstractOptimizer) = string(typeof(opt).name.name)

test_functions = BBOB.list_functions()

b = BBOB.benchmark(
    NelderMead(), BBOB.sphere, [100, 500, 1000], 
)
@test length(b.success_count) == 3

b = BBOB.benchmark(
    NelderMead(), test_functions[1], [100, 500, 1000], Ntrials=20, dimension=2, Δf=1e-6,
)

b = BBOB.benchmark(
    NelderMead(), test_functions[1:3], 100:100:2000, Ntrials=30,
)

b2 = BBOB.benchmark(
    ParticleSwarm(), test_functions[1:3], 100:200:5000, Ntrials=20,
)

plot(b; title = "Benchmark", label = "NelderMead", legend = :outerright)
plot!(b2; label = "ParticleSwarm")

BBOB.compute_CI!(b, 0.1)
