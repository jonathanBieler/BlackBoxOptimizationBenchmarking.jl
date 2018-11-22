using Optim, BlackBoxOptim, NLopt, PyCall

import BlackBoxOptimizationBenchmarking: minimizer, minimum, optimize 
import Base.string

box(D) = fill((-5.5, 5.5),D)
pinit(D) = 10*rand(D).-5

# Define optimize, minimum and minimizer for each optimizer

## NLopt

    mutable struct NLoptOptimMethod 
        s::Symbol
    end
    string(opt::NLoptOptimMethod) = string("NLopt.",opt.s)

    function optimize(NLmeth::NLoptOptimMethod,f,D,run_length)  
        opt = Opt(NLmeth.s, D)
        min_objective!(opt, (p,g) -> f(p) ) #NLopt expect gradient
        maxeval!(opt,run_length)
        xtol_abs!(opt,1e-12)
        xtol_rel!(opt,1e-12)
        ftol_abs!(opt,1e-12)
        ftol_rel!(opt,1e-12)
        lower_bounds!(opt, -5.5*ones(D))
        upper_bounds!(opt, +5.5*ones(D))
        minf,minx,ret = NLopt.optimize(opt, pinit(D))
        return NLmeth, minx, minf
    end
    minimum(mfit::Tuple{NLoptOptimMethod,Vector{Float64},Float64}) = mfit[3]
    minimizer(mfit::Tuple{NLoptOptimMethod,Vector{Float64},Float64}) = mfit[2]

## chain

    mutable struct Chain{T,K}
        first::T
        second::K
        p::Float64
    end
    
    function optimize(m::Chain,f,D,run_length) 
        rl1 = round(Int,m.p*run_length)
        rl2 = run_length - rl1
        
        mfit = optimize(m.first,f,D,run_length)
        xinit = minimizer(mfit)
        mfit = optimize(m.second,f,D,run_length,xinit) 
    end
    
    string(m::Chain) = string(string(m.first)," → ", string(m.second))
    

## python cma

#try

    @pyimport cma
    
    struct PyCMA
    end

    function optimize(m::PyCMA,f,D,run_length)
        es = cma.CMAEvolutionStrategy(pinit(D), 3, Dict("verb_log"=>0,"verb_disp"=>0,"maxfevals"=>run_length))
        mfit = es[:optimize](f)[:result]
        (m,mfit[1],mfit[2])
    end

    string(m::PyCMA) = "PyCMA"
    minimum(mfit::Tuple{PyCMA,Vector{Float64},Float64}) = mfit[3]
    minimizer(mfit::Tuple{PyCMA,Vector{Float64},Float64}) = mfit[2]
  
## scipy

    @pyimport scipy.optimize as scipy_opt

    struct PyMinimize
        method::String
    end

    optimize(m::PyMinimize,f,D,run_length) = (m,scipy_opt.minimize(
        f, pinit(D),method=m.method,
        options = Dict(
            "maxfev"=>run_length,"xatol"=>1e-8,"fatol"=>1e-8,
            "maxiter"=>run_length,"gtol"=>1e-12,
        )
    ))
    minimum(mfit::Tuple{PyMinimize,Dict{Any,Any}}) = mfit[2]["fun"]
    minimizer(mfit::Tuple{PyMinimize,Dict{Any,Any}}) = mfit[2]["x"]
    
    string(m::PyMinimize) = string("Py.",m.method)
    
#end
    
## Optim

    optimize(opt::Optim.AbstractOptimizer,f,D,run_length) =
        Optim.optimize(f, pinit(D), opt, Optim.Options(f_calls_limit=run_length,g_tol=1e-120,iterations=run_length))
                
    optimize(opt::Optim.AbstractOptimizer,f,D,run_length,xinit) =
        Optim.optimize(f, xinit, opt, Optim.Options(f_calls_limit=run_length,g_tol=1e-120,iterations=run_length))
        
    optimize(opt::Optim.SAMIN,f,D,run_length) =
        Optim.optimize(f, fill(-5.5,D), fill(5.5,D), pinit(D), opt, Optim.Options(f_calls_limit=run_length,g_tol=1e-120,iterations=run_length))

    string(opt::Optim.AbstractOptimizer) = string(typeof(opt).name)

    # Optim with restart
    try 
        mutable struct OptimRestart{T}
            opt::T
        end
    catch
    end

    function optimize(opt::OptimRestart,f,D,run_length) 
        fits = [optimize(opt.opt,f,D,round(Int,run_length/20)) for i=1:20]
        mins = [minimum(fit) for fit in fits]
        fits[argmin(mins)]
    end

    string(opt::OptimRestart) = string("Restart-",string(opt.opt))

## BlackBoxOptim

    mutable struct BlackBoxOptimMethod 
        s::Symbol
    end
    string(opt::BlackBoxOptimMethod) = string("BBO.",opt.s)
    
    optimize(method::BlackBoxOptimMethod,f,D,run_length) =
        bboptimize(f; SearchRange=box(D), NumDimensions=D, Method=method.s, MaxFuncEvals=run_length, TraceMode=:silent)

    minimum(mfit::BlackBoxOptim.OptimizationResults) = best_fitness(mfit)
    minimizer(mfit::BlackBoxOptim.OptimizationResults) = best_candidate(mfit)
