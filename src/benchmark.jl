##
import Optim: minimizer, optimize, minimum

mutable struct OptFun
    opt
    fun_idx::Int
end

string(opt::OptFun) = string(string(opt.opt),"\t", enumerate(BBOBFunction)[opt.fun_idx])

function benchmark(op::OptFun, run_lengths, Ntrials, dimensions, Δf)

    f = enumerate(BBOBFunction)[op.fun_idx]
    optimizer = op.opt
    @info("$(string(optimizer))\t $f")

    reached_minium = zeros(Bool,Ntrials,length(run_lengths),length(dimensions))
    distance_to_xopt = zeros(Ntrials,length(run_lengths),length(dimensions))
    fmin = zeros(Ntrials,length(run_lengths),length(dimensions))
    
    t = 0.0
    for j in 1:length(run_lengths), i in 1:Ntrials, k in 1:length(dimensions)
        try
            D = dimensions[k]
            t += @elapsed mfit = optimize(optimizer,f.f,D,run_lengths[j])
            
            reached_minium[i,j,k] = minimum(mfit) < Δf + f.f_opt
            fmin[i,j,k] = minimum(mfit) - f.f_opt 
            distance_to_xopt[i,j,k] =  √sum(abs2.(minimizer(mfit) - f.x_opt[1:D]))
            
        catch err
            @warn(err)
            
            reached_minium[i,j,k] = false
            fmin[i,j,k] = NaN
            distance_to_xopt[i,j,k] = NaN
            @warn(string(optimizer, " failed :", err))
        end

    end
    t /= Ntrials*length(run_lengths)*length(dimensions)
    
    mean(reached_minium,dims=1), mean(distance_to_xopt,dims=1), mean(fmin,dims=1), t
end

benchmark(optimizer::Any, f::BBOBFunction, run_lengths, Ntrials, dimensions, Δf) = benchmark(
    OptFun(optimizer,findfirst(enumerate(BBOBFunction) .== f)), run_lengths, Ntrials, dimensions, Δf
)

benchmark(optimizer::Any, fun_idx::Int, run_lengths, Ntrials, dimensions, Δf) = benchmark(
    OptFun(optimizer,fun_idx), run_lengths, Ntrials, dimensions, Δf
)

#

function benchmark(optimizers::Vector{T}, funcs, run_lengths, Ntrials, dimensions, Δf) where T
    
    mean_succ = zeros(length(optimizers),length(funcs),length(run_lengths),length(dimensions))
    mean_dist, mean_fmin = similar(mean_succ), similar(mean_succ)
    runtime = zeros(length(optimizers),length(funcs))
    
    ops = [OptFun(optimizers[i],funcs[j]) for i=1:length(optimizers), j=1:length(funcs)]
    
    out = pmap(op->benchmark(op, run_lengths, Ntrials, dimensions, Δf), ops)
    
    for i=1:length(optimizers), j=1:length(funcs)
        mean_succ[i,j,:,:], mean_dist[i,j,:,:], mean_fmin[i,j,:,:], runtime[i,j] = out[i,j]
    end
    mean_succ, mean_dist, mean_fmin, runtime
end




##



