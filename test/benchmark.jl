using BBOBFunctions
using Optim, Gadfly, BlackBoxOptim

    
##

import Base.minimum

#Optim
fit(opt::Optim.Optimizer,f,D,run_length) =
    optimize(f, pinit(D), NelderMead(), Optim.Options(iterations=run_length,g_tol=1e-12))
    
minimum(mfit::Optim.OptimizationResults) = mfit.minimum
minimizer(mfit::Optim.OptimizationResults) = mfit.minimizer

# BlackBoxOptim

type BlackBoxOptimMethod 
    s::Symbol
end

box(D) = fill((-5.0, 5.0),D)
pinit(D) = 10*rand(D)-5

fit(method::BlackBoxOptimMethod,f,D,run_length) =
    bboptimize(f; SearchRange=box(D), NumDimensions=D, Method=method.s, MaxFuncEvals=run_length, TraceMode=:silent)

mfit = fit(BlackBoxOptimMethod(:de_rand_1_bin),f.f,D,run_lengths[end])
minimum(mfit) - f.f_opt


minimum(mfit::BlackBoxOptim.OptimizationResults) = best_fitness(mfit)
minimizer(mfit::BlackBoxOptim.OptimizationResults) = best_candidate(mfit)

function run_optimizer(optimizer, run_lengths, Ntrials, f::BBOBFunctions.BBOBFunction, D, Δf)

    reached_minium = zeros(Bool,Ntrials,length(run_lengths))
    distance_to_xopt = zeros(Ntrials,length(run_lengths))
    
    for j in 1:length(run_lengths), i in 1:Ntrials
        mfit = fit(optimizer,f.f,D,run_lengths[j]) 
        reached_minium[i,j] = minimum(mfit) < Δf + f.f_opt
        distance_to_xopt[i,j] =  √sum(abs2(minimizer(mfit) - f.x_opt[1:D]))
    end
    
    mean(reached_minium,1), mean(distance_to_xopt,1)
end

##

D = 2
Ntrials = 20
Δf = 1e-6

run_lengths = round(Int,linspace(2,5000,20))

f = BBOBFunctions.F4

optimizers = [
    NelderMead(), GradientDescent(), BFGS(), 
    BlackBoxOptimMethod(:adaptive_de_rand_1_bin),
    BlackBoxOptimMethod(:xnes),
    BlackBoxOptimMethod(:generating_set_search),
    BlackBoxOptimMethod(:de_rand_2_bin),
]

mean_succ, mean_dist = [], []

for optimizer in optimizers
    m,d = run_optimizer(optimizer, run_lengths, Ntrials, f, D, Δf)
    push!(mean_succ, m)
    push!(mean_dist, d)
end

cols = Colors.distinguishable_colors(length(mean_succ),colorant"red")

plot(
    [layer(x=collect(run_lengths)/D,y=mean_succ[i],Geom.line,Theme(default_color=cols[i])) for i in 1:length(mean_succ)]...,
    Guide.title(f.name),
)

## distance to minimum

plot(
    [layer(x=collect(run_lengths)/D,y=mean_dist[i],Geom.line,Theme(default_color=cols[i])) for i in 1:length(mean_succ)]...,
    Guide.title(f.name),
)

##






