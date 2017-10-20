using Gadfly

#addprocs(2)

## Setup

@everywhere begin
    include(joinpath(Pkg.dir(),"BBOBFunctions","src","benchmark.jl"))
    
    D = 3
    Ntrials = 15
    Δf = 1e-6
    run_lengths = round(Int,linspace(20,20_000,10))
    funcs = 1:14#enumerate(BBOBFunctions.BBOBFunction)
    
    function runopt(op) 
        println(string(op))
        run_optimizer(op, run_lengths, Ntrials, D, Δf)
    end
end

optimizers = [
    OptimRestart(GradientDescent()),
    OptimRestart(NelderMead()),
    CMAESoptim,
    NelderMead(), GradientDescent(), BFGS(), 
    BlackBoxOptimMethod(:adaptive_de_rand_1_bin_radiuslimited),
    BlackBoxOptimMethod(:adaptive_de_rand_1_bin),
    BlackBoxOptimMethod(:xnes),
    BlackBoxOptimMethod(:generating_set_search),
    BlackBoxOptimMethod(:de_rand_2_bin),
]

opt_strings = map(string,optimizers)

## run benchmark

mean_succ = zeros(length(optimizers),length(funcs),length(run_lengths))
mean_dist, mean_fmin = similar(mean_succ), similar(mean_succ)

ops = [OptFun(optimizers[i],funcs[j]) for i=1:length(optimizers), j=1:length(funcs) ]
out = pmap(runopt, ops)

for i=1:length(optimizers), j=1:length(funcs)
    mean_succ[i,j,:], mean_dist[i,j,:], mean_fmin[i,j,:] = out[i,j]
end

## make plots

outdir = joinpath(Pkg.dir(),"BBOBFunctions","data","plots")

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
