using Gadfly
using Distributed
##

run_bench = true

#addprocs(2)

## Setup

@everywhere begin

    using BlackBoxOptimizationBenchmarking
    const BBOB = BlackBoxOptimizationBenchmarking
    include(joinpath(Pkg.dir(),"BlackBoxOptimizationBenchmarking/scripts/optimizers_interface.jl"))

    dimensions = [3 6 16]
    Ntrials = 10
    Δf = 1e-6
    run_lengths = round.(Int,linspace(20,30_000,15))
    #run_lengths = round.(Int,linspace(20,20_000,15))
    funcs = 1:length(enumerate(BBOBFunction))
    
end

optimizers = [
    OptimRestart(NelderMead()),
    NelderMead(), 
    SimulatedAnnealing(),
    BlackBoxOptimMethod(:adaptive_de_rand_1_bin_radiuslimited),
    BlackBoxOptimMethod(:adaptive_de_rand_1_bin),
    Chain(NLoptOptimMethod(:GN_ISRES),NelderMead(),0.9),
    Chain(NLoptOptimMethod(:GN_ESCH), NelderMead(),0.9),
    Chain(NLoptOptimMethod(:GD_STOGO),NelderMead(),0.9),
    BlackBoxOptimMethod(:xnes),
    BlackBoxOptimMethod(:generating_set_search),
    BlackBoxOptimMethod(:de_rand_2_bin),
    #BlackBoxOptimMethod(:resampling_memetic_search),#poor performance 
    PyMinimize("Nelder-Mead"),
    PyCMA(),
]

opt_strings = map(string,optimizers)

## run benchmark

if run_bench

    mean_succ, mean_dist, mean_fmin, runtime = BBOB.benchmark(
        optimizers, funcs, run_lengths, Ntrials, dimensions, Δf,
    )

    @info("Done, saving data and making plots.")

    # save output

    outdir = joinpath(Pkg.dir(),"BlackBoxOptimizationBenchmarking","data")

    writedlm(joinpath(outdir,"mean_succ.txt"),mean_succ)
    writedlm(joinpath(outdir,"mean_dist.txt"),mean_dist)
    writedlm(joinpath(outdir,"mean_fmin.txt"),mean_fmin)
    writedlm(joinpath(outdir,"runtime.txt"),runtime)

else

    # load output
    outdir = joinpath(Pkg.dir(),"BlackBoxOptimizationBenchmarking","data")

    s = (length(optimizers),length(funcs),length(run_lengths),length(dimensions))
    mean_succ = reshape( readdlm(joinpath(outdir,"mean_succ.txt")), s) 
    mean_dist = reshape( readdlm(joinpath(outdir,"mean_dist.txt")), s)
    mean_fmin = reshape( readdlm(joinpath(outdir,"mean_fmin.txt")), s)
    runtime = readdlm(joinpath(outdir,"runtime.txt"))

end

## make plots

outdir = joinpath(Pkg.dir(),"BlackBoxOptimizationBenchmarking","data","plots")
cols = Colors.distinguishable_colors(size(mean_succ,1)+1,colorant"white")[2:end]

#all together

idx = sortperm([mean(mean_succ[i,:,end,:]) for i=1:length(optimizers)],rev=true)

p = plot(
    [layer(
        x=collect(run_lengths),y=mean(mean_succ[i,:,:,:],(1,3)),Geom.line,
        Theme(default_color=cols[i], line_style= mod(i,2)==0 ? :solid : :dash, line_width=2pt)
    ) for i in 1:size(mean_succ,1)]...,
    Coord.cartesian(ymax=1),
    Guide.title("All functions"), Guide.xlabel("Run Length"), Guide.ylabel("Success rate"),
    Guide.manual_color_key("", opt_strings[idx], cols[idx])
)

#

draw(PDF(joinpath(outdir,"mean_succ.pdf"),22cm,14cm),p)
draw(PNG(joinpath(outdir,"mean_succ.png"),dpi=150,16cm,10cm),p)

# each dimension
for (k,D) in zip(1:length(dimensions), dimensions)

    idx = sortperm([mean(mean_succ[i,:,end,k]) for i=1:length(optimizers)],rev=true)

    p = plot(
        [layer(
            x=collect(run_lengths),y=mean(mean_succ[i,:,:,k],1),Geom.line,
            Theme(default_color=cols[i], line_style= mod(i,2)==0 ? :solid : :dash, line_width=2pt)
        ) for i in 1:size(mean_succ,1)]...,
        Coord.cartesian(ymax=1),
        Guide.title("All functions, D: $D"), Guide.xlabel("Run Length"), Guide.ylabel("Success rate"), 
        Guide.manual_color_key("", opt_strings[idx], cols[idx])
    )
    draw(PDF(joinpath(outdir,"per_dimension","mean_succ_$(D).pdf"),22cm,14cm),p)
    draw(PNG(joinpath(outdir,"per_dimension","mean_succ_$(D).png"),dpi=150,16cm,10cm),p)

end

## f min

idx = sortperm([exp(mean(log.(mean_fmin[i,:,end,:]+eps()))) for i=1:length(optimizers)],rev=true)

# all together
p = plot(
    [layer(
        x=collect(run_lengths),y=exp.(mean(log.(mean_fmin[i,:,:,:]+eps()),(1,3))),Geom.line,
        Theme(default_color=cols[i],line_width=2pt)
    ) 
    for i in 1:size(mean_succ,1)]...,
    layer(yintercept=[Δf],Geom.hline,Theme(default_color=colorant"gray")),
    Guide.title("All functions"), Guide.xlabel("Run Length"), Guide.ylabel("fmin"), 
    Guide.manual_color_key("", opt_strings[idx], cols[idx]),
#    Scale.x_log10,
    Scale.y_log10,
)

#
draw(PDF(joinpath(outdir,"mean_fmin.pdf"),22cm,14cm),p)
draw(PNG(joinpath(outdir,"mean_fmin.png"),dpi=150,16cm,10cm),p)

for (k,D) in zip(1:length(dimensions), dimensions)

    idx = sortperm([exp.(mean(log.(mean_fmin[i,:,end,k]+eps()))) for i=1:length(optimizers)],rev=true)

    p = plot(
        [layer(
            x=collect(run_lengths)/D,y=exp(mean(log.(mean_fmin[i,:,:,k]+eps()),1)),Geom.line,
            Theme(default_color=cols[i],line_width=2pt)
        ) for i in 1:size(mean_succ,1)]...,
        layer(yintercept=[Δf],Geom.hline,Theme(default_color=colorant"gray")),
        Guide.title("All functions"), Guide.xlabel("Run Length / D"), Guide.ylabel("fmin"), 
        Guide.manual_color_key("", opt_strings[idx], cols[idx]),
    #    Scale.x_log10,
        Scale.y_log10,
    )

    draw(PDF(joinpath(outdir,"per_dimension","mean_fmin_$(D).pdf"),20cm,12cm),p)
    draw(PNG(joinpath(outdir,"per_dimension","mean_fmin_$(D).png"),dpi=150,16cm,10cm),p)

end

## distance to minimizer

for (k,D) in zip(1:length(dimensions), dimensions)

    p = plot(
        [layer(x=collect(run_lengths)/D,y=median(mean_dist[i,:,:,k],1),Geom.line,Theme(default_color=cols[i])) for i in 1:size(mean_succ,1)]...,
        Guide.title("All functions"), Guide.xlabel("Run Length / D"), Guide.ylabel("Distance to xmin"), 
        Guide.manual_color_key("", opt_strings, cols),
    #    Scale.x_log10,
        Scale.y_log10,
    )

#    draw(PDF(joinpath(outdir,"per_dimension","mean_dist_$(D).pdf"),20cm,12cm),p)
#    draw(SVG(joinpath(outdir,"per_dimension","mean_dist_$(D).svg"),20cm,12cm),p)

end

## runtime

idx = sortperm(mean(runtime,2)[:],rev=true)

p = plot(
    layer(
        y=opt_strings[idx], x=(mean(runtime,2)/Base.minimum(mean(runtime,2)))[idx],
        color = cols[idx],
        Geom.bar(orientation=:horizontal)
    ),
    Guide.title("All functions"), Guide.ylabel(""), 
    Guide.xlabel("Relative Run Time",orientation=:horizontal), 
    Scale.x_log10,
#    Scale.y_log10,
)

draw(PDF(joinpath(outdir,"runtime.pdf"),20cm,12cm),p)
draw(PNG(joinpath(outdir,"runtime.png"),dpi=150,16cm,10cm),p)
p

##
