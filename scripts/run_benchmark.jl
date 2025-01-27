##

using BlackBoxOptimizationBenchmarking, Plots, Optimization, Memoize, Statistics
import BlackBoxOptimizationBenchmarking.Chain
const BBOB = BlackBoxOptimizationBenchmarking

using OptimizationBBO, OptimizationOptimJL, OptimizationEvolutionary, OptimizationNLopt
using OptimizationMetaheuristics, OptimizationNOMAD

chain = (t; isboxed=false) -> Chain(
    BenchmarkSetup(t, isboxed = isboxed),
    BenchmarkSetup(NelderMead(), isboxed = false),
    0.9
)

test_functions = BBOB.list_functions()

@memoize run_bench(algo) = BBOB.benchmark(setup[algo], test_functions, run_length, Ntrials=30, dimension = 3)

## test one

setup = Dict(
    "NelderMead" => NelderMead(),
    #Optim.BFGS(),
    "NLopt.GN_MLSL_LDS" => chain(NLopt.GN_MLSL_LDS(), isboxed=true),
    "NLopt.GN_CRS2_LM()" => chain(NLopt.GN_CRS2_LM(), isboxed=true),
    "NLopt.GN_DIRECT()" => chain(NLopt.GN_DIRECT(), isboxed=true),
    "NLopt.GN_ESCH()"  => chain(NLopt.GN_ESCH(), isboxed=true),
    "OptimizationEvolutionary.GA()" => chain(OptimizationEvolutionary.GA()),
    "OptimizationEvolutionary.DE()" => chain(OptimizationEvolutionary.DE()),
    "OptimizationEvolutionary.ES()" => chain(OptimizationEvolutionary.ES()),
    "Optim.SAMIN" => chain(SAMIN(verbosity=0), isboxed=true),
    "BBO_adaptive_de_rand_1_bin" => chain(BBO_adaptive_de_rand_1_bin(), isboxed=true),
    "BBO_separable_nes" => chain(BBO_separable_nes(), isboxed=true),
    "BBO_xnes" => chain(BBO_xnes(), isboxed=true),
    # "NOMADOpt" => chain(NOMADOpt()), too much printing
    "OptimizationMetaheuristics.ECA" => chain(OptimizationMetaheuristics.ECA(), isboxed=true),
    "OptimizationMetaheuristics.CGSA" => chain(OptimizationMetaheuristics.CGSA(), isboxed=true),
    "OptimizationMetaheuristics.DE" => chain(OptimizationMetaheuristics.DE(), isboxed=true),
    #chain(BBO_adaptive_de_rand_1_bin_radiuslimited(), isboxed=true), # same as BBO_adaptive_de_rand_1_bin
)

@time b = BBOB.benchmark(
    setup["OptimizationMetaheuristics.CGSA"], test_functions[1:3], 100:500:5000, Ntrials=10, dimension = 3
)

plot(b)

## try all

results = Array{BBOB.BenchmarkResults}(undef, length(setup))

run_length = round.(Int, 10 .^ LinRange(1,4,25))

Threads.@threads for (i,algo) in collect(enumerate(keys(setup)))
    results[i] = run_bench(algo)
end

## plot

labels = collect(keys(setup))
idx = sortperm([b.success_rate[end] for b in results], rev=true)

p = plot(xscale = :log10, legend = :outerright, size = (800,600), margin=10Plots.px)
for i in idx
    plot!(results[i], label = labels[i], showribbon=false, lw=2)
end
p

## make heatmap per function

success_rate_per_function = reduce(hcat, b.success_rate_per_function for b in results)

idx = sortperm(mean(success_rate_per_function, dims=1)[:], rev=false)
idxfunc = sortperm(mean(success_rate_per_function, dims=2)[:], rev=true)
idxfunc = 1:length(test_functions)

heatmap(
    string.(test_functions)[idxfunc], labels[idx], success_rate_per_function[idxfunc, idx]',
    cmap = :RdYlGn,
    xticks = :all,
    xrotation = 45,
)

##

bar(getfield.(results, :runtime))

## run a single problem and plot solution

Δf = 1e-5
f = test_functions[3]
algo = "NLopt.GN_CRS2_LM()"
sol = [BBOB.solve_problem(setup[algo], f, 3, 5_000) for in in 1:10]

@info [sol.objective < Δf + f.f_opt for sol in sol]

p = plot(f, size = (600,600), zoom = 1)
for sol in sol
    scatter!(sol.u[1:1], sol.u[2:2], label="", c="blue", marker = :xcross)
end
p

##
