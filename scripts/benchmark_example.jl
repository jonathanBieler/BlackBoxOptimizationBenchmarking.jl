using BlackBoxOptimizationBenchmarking, Optim
import BlackBoxOptimizationBenchmarking: minimizer, minimum, optimize 
import Base.string

const BBOB = BlackBoxOptimizationBenchmarking

box(D) = fill((-5.5, 5.5),D)
pinit(D) = 10*rand(D).-5

optimize(opt::Optim.SAMIN,f,D,run_length) =
    Optim.optimize(f, fill(-5.5,D), fill(5.5,D), pinit(D), opt, Optim.Options(f_calls_limit=run_length,g_tol=1e-120,iterations=run_length))

run_length = 30_000
dimensions = 3
Ntrials = 20
Δf = 1e-6

# test on Sphere
f = BlackBoxOptimizationBenchmarking.F1
success_rate, x_dist, f_dist, runtime = BBOB.benchmark(Optim.SAMIN(), f, run_length, Ntrials, dimensions, Δf)

#test on all functions
perf = [(f=>BBOB.benchmark(Optim.SAMIN(), f, run_length, Ntrials, dimensions, Δf)[1][1]) for f in enumerate(BBOBFunction)]
    