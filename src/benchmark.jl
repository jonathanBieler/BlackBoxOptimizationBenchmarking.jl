using BBOBFunctions
using Optim, BlackBoxOptim
import Base: string
    
##

include("optimizers_interface.jl")

##

type OptFun
    opt
    fun_idx::Int
end

string(opt::OptFun) = string(string(opt.opt),"\t", enumerate(BBOBFunctions.BBOBFunction)[opt.fun_idx])

function run_optimizer(op::OptFun, run_lengths, Ntrials, D, Δf)

    f = enumerate(BBOBFunctions.BBOBFunction)[op.fun_idx]
    optimizer = op.opt

    reached_minium = zeros(Bool,Ntrials,length(run_lengths))
    distance_to_xopt = zeros(Ntrials,length(run_lengths))
    fmin = zeros(Ntrials,length(run_lengths))
    
    for j in 1:length(run_lengths), i in 1:Ntrials
        mfit = fit(optimizer,f.f,D,run_lengths[j]) 
        reached_minium[i,j] = minimum(mfit) < Δf + f.f_opt
        fmin[i,j] = minimum(mfit) - f.f_opt 
        distance_to_xopt[i,j] =  √sum(abs2.(minimizer(mfit) - f.x_opt[1:D]))
    end
    
    mean(reached_minium,1), mean(distance_to_xopt,1), mean(fmin,1)
end

##



