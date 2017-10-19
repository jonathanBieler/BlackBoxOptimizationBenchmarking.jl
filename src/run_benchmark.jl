using Gadfly

#addprocs(2)

## Setup

@everywhere begin
    include(joinpath(Pkg.dir(),"BBOBFunctions","src","benchmark.jl"))
    
    D = 3
    Ntrials = 10
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

## run benchmark

mean_succ = zeros(length(optimizers),length(funcs),length(run_lengths))
mean_dist, mean_fmin = similar(mean_succ), similar(mean_succ)

for fi = 1:length(funcs)
    warn(enumerate(BBOBFunctions.BBOBFunction)[fi])

    ops = [OptFun(optimizers[i],funcs[fi]) for i=1:length(optimizers)]
    out = pmap(runopt, ops)
#    local_out = map(run, ops)
    
    for oi in 1:length(optimizers)
        mean_succ[oi,fi,:], mean_dist[oi,fi,:], mean_fmin[oi,fi,:] = out[oi]
    end
end

## make plots

outdir = joinpath(Pkg.dir(),"BBOBFunctions","data","plots")

opt_strings = ["GradientDescent-R","CMAES","Simplex","GD","BFGS","a_de_rand_1_bin","xnes","generating_set_search","de_rand_2_bin"]
cols = Colors.distinguishable_colors(size(mean_succ,1),colorant"red")

p = plot(
    [layer(x=collect(run_lengths)/D,y=mean(mean_succ[i,:,:],1),Geom.line,Theme(default_color=cols[i])) for i in 1:size(mean_succ,1)]...,
    Guide.title("All functions"),
    Guide.manual_color_key("", opt_strings, cols)
)
draw(PDF(joinpath(outdir,"mean_succ.pdf"),16cm,12cm),p)

## f min

p = plot(
    [layer(x=collect(run_lengths)/D,y=median(mean_fmin[i,:,:],1),Geom.line,Theme(default_color=cols[i])) for i in 1:size(mean_succ,1)]...,
    layer(yintercept=[Δf],Geom.hline,Theme(default_color=colorant"gray")),
    Guide.title("All functions"),
    Guide.manual_color_key("", opt_strings, cols),
#    Scale.x_log10,
    Scale.y_log10,
)

draw(PDF(joinpath(outdir,"mean_fmin.pdf"),16cm,12cm),p)

## distance to minimizer

p = plot(
    [layer(x=collect(run_lengths)/D,y=median(mean_dist[i,:,:],1),Geom.line,Theme(default_color=cols[i])) for i in 1:size(mean_succ,1)]...,
    Guide.title("All functions"),
    Guide.manual_color_key("", opt_strings, cols),
#    Scale.x_log10,
    Scale.y_log10,
)

draw(PDF(joinpath(outdir,"mean_dist.pdf"),16cm,12cm),p)

##
