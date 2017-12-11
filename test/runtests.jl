using Base.Test

# benchmark uses several workers if available
@everywhere begin

    using BlackBoxOptimizationBenchmarking, Optim
    import BlackBoxOptimizationBenchmarking: minimizer, minimum, optimize
    const BBOB = BlackBoxOptimizationBenchmarking
    
    #interface for Optim
    pinit(D) = 10*rand(D)-5
    optimize(opt::Optim.Optimizer,f,D,run_length) =
        Optim.optimize(f, pinit(D), opt, Optim.Options(f_calls_limit=run_length,g_tol=1e-12))
        
    minimum(mfit::Optim.OptimizationResults) = mfit.minimum
    minimizer(mfit::Optim.OptimizationResults) = mfit.minimizer

    string(opt::Optim.Optimizer) = string(typeof(opt).name)
end

b = BBOB.benchmark(
    NelderMead(), BBOB.F1, [100,500,1000], 20, 2, 1e-6,
)

b = BBOB.benchmark(
    NelderMead(), 1, [100,500,1000], 20, 2, 1e-6,
)

b = BBOB.benchmark(
    [NelderMead(), GradientDescent()], [1, 2], [100,500,1000], 20, 2, 1e-6,
)

