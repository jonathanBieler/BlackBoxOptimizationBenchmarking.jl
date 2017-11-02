using Gadfly

#addprocs(2)

## Setup

@everywhere begin
    include(joinpath(Pkg.dir(),"BBOBFunctions","src","benchmark.jl"))
    
    D = 5
    Ntrials = 30
    Δf = 1e-6
    run_lengths = round.(Int,linspace(20,40_000,30))
    funcs = 1:length(enumerate(BBOBFunctions.BBOBFunction))
    
    function runopt(op) 
        println(string(op))
        run_optimizer(op, run_lengths, Ntrials, D, Δf)
    end
end

optimizers = [
    OptimRestart(GradientDescent()),
    OptimRestart(NelderMead()),
    CMAESoptim,
    NelderMead(), GradientDescent(), 
    BlackBoxOptimMethod(:adaptive_de_rand_1_bin_radiuslimited),
    BlackBoxOptimMethod(:adaptive_de_rand_1_bin),
    NLoptOptimMethod(:LN_BOBYQA),
    #NLoptOptimMethod(:LN_PRAXIS), #this guy often fails
    #NLoptOptimMethod(:LN_SBPLX), # this guy segfault
    BlackBoxOptimMethod(:xnes),
    BlackBoxOptimMethod(:generating_set_search),
    BlackBoxOptimMethod(:de_rand_2_bin),
    BlackBoxOptimMethod(:resampling_memetic_search),
]

opt_strings = map(string,optimizers)

## run benchmark

mean_succ = zeros(length(optimizers),length(funcs),length(run_lengths))
mean_dist, mean_fmin = similar(mean_succ), similar(mean_succ)

runtime = zeros(length(optimizers),length(funcs))

ops = [OptFun(optimizers[i],funcs[j]) for i=1:length(optimizers), j=1:length(funcs) ]
out = pmap(runopt, ops)

for i=1:length(optimizers), j=1:length(funcs)
    mean_succ[i,j,:], mean_dist[i,j,:], mean_fmin[i,j,:], runtime[i,j] = out[i,j]
end

## save output
outdir = joinpath(Pkg.dir(),"BBOBFunctions","data")

writedlm(joinpath(outdir,"mean_succ.txt"),mean_succ)
writedlm(joinpath(outdir,"mean_dist.txt"),mean_dist)
writedlm(joinpath(outdir,"mean_fmin.txt"),mean_fmin)

## make plots

outdir = joinpath(Pkg.dir(),"BBOBFunctions","data","plots")
cols = Colors.distinguishable_colors(size(mean_succ,1),colorant"red")

p = plot(
    [layer(x=collect(run_lengths)/D,y=mean(mean_succ[i,:,:],1),Geom.line,Theme(default_color=cols[i])) for i in 1:size(mean_succ,1)]...,
    Guide.title("All functions"), Guide.xlabel("Run Length / D"), Guide.ylabel("Success rate"), 
    Guide.manual_color_key("", opt_strings, cols)
)
draw(PDF(joinpath(outdir,"mean_succ.pdf"),18cm,12cm),p)
draw(SVG(joinpath(outdir,"mean_succ.svg"),18cm,12cm),p)

## f min

p = plot(
    [layer(x=collect(run_lengths)/D,y=median(mean_fmin[i,:,:],1),Geom.line,Theme(default_color=cols[i])) for i in 1:size(mean_succ,1)]...,
    layer(yintercept=[Δf],Geom.hline,Theme(default_color=colorant"gray")),
    Guide.title("All functions"), Guide.xlabel("Run Length / D"), Guide.ylabel("fmin"), 
    Guide.manual_color_key("", opt_strings, cols),
#    Scale.x_log10,
    Scale.y_log10,
)

draw(PDF(joinpath(outdir,"mean_fmin.pdf"),18cm,12cm),p)
draw(SVG(joinpath(outdir,"mean_fmin.svg"),18cm,12cm),p)

## distance to minimizer

p = plot(
    [layer(x=collect(run_lengths)/D,y=median(mean_dist[i,:,:],1),Geom.line,Theme(default_color=cols[i])) for i in 1:size(mean_succ,1)]...,
    Guide.title("All functions"), Guide.xlabel("Run Length / D"), Guide.ylabel("Distance to xmin"), 
    Guide.manual_color_key("", opt_strings, cols),
#    Scale.x_log10,
    Scale.y_log10,
)

draw(PDF(joinpath(outdir,"mean_dist.pdf"),18cm,12cm),p)
draw(SVG(joinpath(outdir,"mean_dist.svg"),18cm,12cm),p)

## runtime

p = plot(
    layer(x=opt_strings, y=mean(runtime,2)/minimum(mean(runtime,2)),Geom.bar),
    Guide.title("All functions"), Guide.xlabel(""), Guide.ylabel("Run Time"), 
    #Scale.x_log10,
    #Scale.y_log10,
)

draw(PDF(joinpath(outdir,"runtime.pdf"),18cm,12cm),p)
draw(SVG(joinpath(outdir,"runtime.svg"),18cm,12cm),p)

##
