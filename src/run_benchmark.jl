using Gadfly


#addprocs(2)

@everywhere begin
    
    include("benchmark.jl")
    
    D = 3
    Ntrials = 20
    Δf = 1e-6
    run_lengths = round(Int,linspace(20,10_000,10))
    funcs = 1:14#enumerate(BBOBFunctions.BBOBFunction)
    
    runopt(op) = run_optimizer(op, run_lengths, Ntrials, D, Δf)
end

optimizers = [
    OptimRestart(GradientDescent()),
    CMAESoptim,
    NelderMead(), GradientDescent(), BFGS(), 
    BlackBoxOptimMethod(:adaptive_de_rand_1_bin),
    BlackBoxOptimMethod(:xnes),
    BlackBoxOptimMethod(:generating_set_search),
    BlackBoxOptimMethod(:de_rand_2_bin),
]


##

mean_succ = zeros(length(optimizers),length(funcs),length(run_lengths))
mean_dist, mean_fmin = similar(mean_succ), similar(mean_succ)

for fi = 1:length(funcs)
    warn(funcs[fi])

#    for oi in 1:length(optimizers)
#        mean_succ[oi,fi,:], mean_dist[oi,fi,:], mean_fmin[oi,fi,:] = run_optimizer(optimizers[oi], run_lengths, Ntrials, funcs[fi], D, Δf)
#    end

    ops = [OptFun(optimizers[i],funcs[fi]) for i=1:length(optimizers)]

    
    out = pmap(runopt, ops)
#    local_out = map(run, ops)
    
    for oi in 1:length(optimizers)
        mean_succ[oi,fi,:], mean_dist[oi,fi,:], mean_fmin[oi,fi,:] = out[oi]
    end
end


##

# plot

opt_strings = ["GradientDescent-R","CMAES","Simplex","GD","BFGS","a_de_rand_1_bin","xnes","generating_set_search","de_rand_2_bin"]
cols = Colors.distinguishable_colors(size(mean_succ,1),colorant"red")

plot(
    [layer(x=collect(run_lengths)/D,y=mean(mean_succ[i,:,:],1),Geom.line,Theme(default_color=cols[i])) for i in 1:size(mean_succ,1)]...,
    Guide.title("All functions"),
    Guide.manual_color_key("", opt_strings, cols)
)

## f min

plot(
    [layer(x=collect(run_lengths)/D,y=mean(mean_fmin[i,:,:],1),Geom.line,Theme(default_color=cols[i])) for i in 1:size(mean_succ,1)]...,
    layer(yintercept=[Δf],Geom.hline,Theme(default_color=colorant"gray")),
    Guide.title("All functions"),
    Guide.manual_color_key("", opt_strings, cols),
#    Scale.x_log10,
#    Scale.y_log10,
)


## distance to minimizer

plot(
    [layer(x=collect(run_lengths)/D,y=mean(mean_dist[i,:,:],1),Geom.line,Theme(default_color=cols[i])) for i in 1:size(mean_succ,1)]...,
#    layer(yintercept=[Δf],Geom.hline,Theme(default_color=colorant"gray")),
    Guide.title("All functions"),
    Guide.manual_color_key("", opt_strings, cols),
#    Scale.x_log10,
    Scale.y_log10,
)

##
