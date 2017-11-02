import Base: minimum, string

# Define fit, minimum and minimizer for each optimizer

## NLopt

type NLoptOptimMethod 
    s::Symbol
end
string(opt::NLoptOptimMethod) = string("NLopt.",opt.s)

function fit(NLmeth::NLoptOptimMethod,f,D,run_length)  
    opt = Opt(NLmeth.s, D)
    min_objective!(opt, (p,g) -> f(p) ) #NLopt expect gradient
    maxeval!(opt,run_length)
    minf,minx,ret = NLopt.optimize(opt, pinit(D))
    return NLmeth, minx, minf
end
minimum(mfit::Tuple{NLoptOptimMethod,Array{Float64,1},Float64}) = mfit[3]
minimizer(mfit::Tuple{NLoptOptimMethod,Array{Float64,1},Float64}) = mfit[2]

#opt = Opt(:LN_BOBYQA, 3)
#min_objective!(opt,(p,g)->BBOBFunctions.F1.f(p))
#maxeval!(opt,1000)
#minf,minx,ret = NLopt.optimize(opt, pinit(D))

## my cmaes

    include("cmaes.jl")
    using CMAES

    type CMAESoptim end

    fit(::Type{CMAESoptim},f,D,run_length) = CMAES.cmaes(f, pinit(D), 3.0, run_length, round(Int, 3 + floor(3log(D))))
    minimum(mfit::Tuple{Array{Float64,1},Float64}) = mfit[2]
    minimizer(mfit::Tuple{Array{Float64,1},Float64}) = mfit[1]

## Optim

    fit(opt::Optim.Optimizer,f,D,run_length) =
        Optim.optimize(f, pinit(D), NelderMead(), Optim.Options(f_calls_limit=run_length,g_tol=1e-12))
        
    minimum(mfit::Optim.OptimizationResults) = mfit.minimum
    minimizer(mfit::Optim.OptimizationResults) = mfit.minimizer

    string(opt::Optim.Optimizer) = string(typeof(opt).name)

    # Optim with restart
    try 
        type OptimRestart{T}
            opt::T
        end
    end

    function fit(opt::OptimRestart,f,D,run_length) 
        fits = [fit(opt.opt,f,D,round(Int,run_length/20)) for i=1:20]
        mins = [minimum(fit) for fit in fits]
        fits[indmin(mins)]
    end

    string(opt::OptimRestart) = string("R-",string(opt.opt))

## BlackBoxOptim

    type BlackBoxOptimMethod 
        s::Symbol
    end
    string(opt::BlackBoxOptimMethod) = string("BBO.",opt.s)
    
    box(D) = fill((-5.0, 5.0),D)
    pinit(D) = 10*rand(D)-5

    fit(method::BlackBoxOptimMethod,f,D,run_length) =
        bboptimize(f; SearchRange=box(D), NumDimensions=D, Method=method.s, MaxFuncEvals=run_length, TraceMode=:silent)

    minimum(mfit::BlackBoxOptim.OptimizationResults) = best_fitness(mfit)
    minimizer(mfit::BlackBoxOptim.OptimizationResults) = best_candidate(mfit)
