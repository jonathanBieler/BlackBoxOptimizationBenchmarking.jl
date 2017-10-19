import Base.minimum

# Define fit, minimum and minimizer for each optimizer

## my cmaes

    include("cmaes.jl")
    using CMAES

    type CMAESoptim end

    fit(::Type{CMAESoptim},f,D,run_length) = CMAES.cmaes(f, pinit(D), 1, run_length, round(Int, 3 + floor(3log(D))))
    minimum(mfit::Tuple{Array{Float64,1},Float64}) = mfit[2]
    minimizer(mfit::Tuple{Array{Float64,1},Float64}) = mfit[1]

## Optim

    fit(opt::Optim.Optimizer,f,D,run_length) =
        optimize(f, pinit(D), NelderMead(), Optim.Options(f_calls_limit=run_length,g_tol=1e-12))
        
    minimum(mfit::Optim.OptimizationResults) = mfit.minimum
    minimizer(mfit::Optim.OptimizationResults) = mfit.minimizer

    # Optim with restart
    try 
        type OptimRestart{T<:Optim.Optimizer}
            opt::T
        end
    end

    function fit(opt::OptimRestart,f,D,run_length) 
        fits = [fit(opt.opt,f,D,round(Int,run_length/20)) for i=1:20]
        mins = [minimum(fit) for fit in fits]
        fits[indmin(mins)]
    end

## BlackBoxOptim

    type BlackBoxOptimMethod 
        s::Symbol
    end
    
    box(D) = fill((-5.0, 5.0),D)
    pinit(D) = 10*rand(D)-5

    fit(method::BlackBoxOptimMethod,f,D,run_length) =
        bboptimize(f; SearchRange=box(D), NumDimensions=D, Method=method.s, MaxFuncEvals=run_length, TraceMode=:silent)

    minimum(mfit::BlackBoxOptim.OptimizationResults) = best_fitness(mfit)
    minimizer(mfit::BlackBoxOptim.OptimizationResults) = best_candidate(mfit)

    