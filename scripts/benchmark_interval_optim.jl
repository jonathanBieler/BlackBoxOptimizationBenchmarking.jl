using Gadfly, Colors, BlackBoxOptimizationBenchmarking
using Distributed, Statistics

bbob_dir = dirname(pathof(BlackBoxOptimizationBenchmarking))

using BlackBoxOptimizationBenchmarking
const BBOB = BlackBoxOptimizationBenchmarking
bbob_dir = dirname(pathof(BBOB))
include(joinpath(bbob_dir,"../scripts/optimizers_interface.jl"))

using IntervalArithmetic, IntervalOptimisation

struct IntervalOptimisationMethod
end
string(opt::IntervalOptimisationMethod) = "IntervalOptimisation"

function optimize(NLmeth::IntervalOptimisationMethod,f,D,run_length)

    minf, minx = minimise(f, IntervalBox(-5.5..5.5, D), f_calls_limit=run_length)
    return NLmeth, Vector(mean(mid.(minx))), mid(minf)
end
minimum(  mfit::Tuple{IntervalOptimisationMethod, Vector{Float64}, Float64}) = mfit[3]
minimizer(mfit::Tuple{IntervalOptimisationMethod, Vector{Float64}, Float64}) = mfit[2]

dimensions = 6
Ntrials = 30
Δf = 1e-6

c = opt -> Chain(opt, NelderMead(), 0.9)

funcs = 1:12 #length(enumerate(BBOBFunction))

run_lengths = round.(Int, range(20, stop=3_000, length=10))

optimizers = [
    c(BlackBoxOptimMethod(:adaptive_de_rand_1_bin)),
    NelderMead(),
    c(IntervalOptimisationMethod()),
]
opt_strings = map(string, optimizers)

#error("")

mean_succ, mean_dist, mean_fmin, runtime = BBOB.benchmark(
    optimizers, funcs, run_lengths, Ntrials, dimensions, Δf,
)

cols = Colors.distinguishable_colors(size(mean_succ,1)+1,colorant"white")[2:end]
line_style = i -> mod(i,2)==0 ? [:solid] : [:dash]

#all together

idx = sortperm([mean(mean_succ[i,:,end,:]) for i=1:length(optimizers)],rev=true)

p = plot(
    [layer(
        x=collect(run_lengths),y=mean(mean_succ[i,:,:,:],dims=(1,3)),Geom.line,
        Theme(default_color=cols[i], line_style=line_style(i), line_width=2pt)
    ) for i in 1:size(mean_succ,1)]...,
    Coord.cartesian(ymax=1),
    Guide.title("All functions"), Guide.xlabel("Run Length"), Guide.ylabel("Success rate"),
    Guide.manual_color_key("", opt_strings[idx], cols[idx]),
    Theme(background_color="white"),
)

draw(PNG("mean_succ.png",16cm,10cm),p)

## runtime

idx = sortperm(mean(runtime,dims=2)[:],rev=true)

p = plot(
    layer(
        y=opt_strings[idx], x=(mean(runtime,dims=2)/Base.minimum(mean(runtime,dims=2)))[idx],
        color = cols[idx],
        Geom.bar(orientation=:horizontal)
    ),
    Guide.title("All functions"), Guide.ylabel(""), 
    Guide.xlabel("Relative Run Time",orientation=:horizontal), 
    Scale.x_log10,
#    Scale.y_log10,
)
