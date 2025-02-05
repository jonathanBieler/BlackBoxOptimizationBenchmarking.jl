##

mutable struct BenchmarkResults
    run_length::Vector{Int}
    success_count::Vector{Int}
    success_rate::Vector{Float64}
    success_rate_qlow::Vector{Float64}
    success_rate_qhigh::Vector{Float64}
    distance_to_minimizer::Vector{Float64}
    minimum::Vector{Float64}
    runtime::Float64
    Neffective::Int
    callcount::Vector{Float64}
    success_rate_per_function::Vector{Float64}
end
BenchmarkResults(; kwargs...) = BenchmarkResults(values(kwargs)...)

show(io::IO, b::BenchmarkResults) =  begin
    println(io, "BenchmarkResults :")
    print(io, "Run length : ")
    show(IOContext(io, :limit => true, :compact => true), b.run_length)
    print(io, "\nSuccess rate : ")
    show(IOContext(io, :limit => true, :compact => true), b.success_rate)
end

# https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Clopper%E2%80%93Pearson_interval
compute_CI(rate::Real, N::Int, q::Real) = quantile(Beta(0.1 + round(Int, N*rate), 0.1 + N - round(Int, N*rate)), q)

function compute_CI!(b::BenchmarkResults, CI_quantile)
    qs = compute_CI(b.success_rate, b.Neffective, CI_quantile)
    b.success_rate_qlow = qs.success_rate_qlow
    b.success_rate_qhigh = qs.success_rate_qhigh
    b
end

function compute_CI(success_rate::Vector{Float64}, Neffective, CI_quantile::Real)
    (;
    success_rate_qlow = compute_CI.(success_rate, Neffective, CI_quantile),
    success_rate_qhigh = compute_CI.(success_rate, Neffective, 1-CI_quantile)
    )
end

pinit(D) = 10*rand(D).-5

mutable struct BenchmarkSetup{T}
    method::T
    isboxed::Bool
end
BenchmarkSetup(method::T; isboxed=false) where T = BenchmarkSetup{T}(method, isboxed)

show(io::IO, b::BenchmarkSetup{T}) where {T} =  begin
    println(io, nameof(T))
end

# FunctionCallsCounter : keep count of how many time our function is called

mutable struct FunctionCallsCounter
    f::Function
    @atomic count::Int
end
FunctionCallsCounter(f::Function) = FunctionCallsCounter(f,0)

function (f::FunctionCallsCounter)(args)
    @atomic f.count += 1
    f.f(args)
end

# Allow to chain two optimizers
mutable struct Chain{T, K}
    first::T
    second::K
    p::Float64
end

show(io::IO, b::Chain{T, K}) where {T,K} =  begin
    T1 = nameof(typeof(b.first.method))
    T2 = nameof(typeof(b.second.method))
    println(io, "Chain($(T1) → $(T2))")
end

function solve_problem(m::Chain, f, D::Int, run_length::Int) 
    rl1 = round(Int, m.p*run_length)
    rl2 = run_length - rl1
    
    sol = solve_problem(m.first, f, D, rl1)
    xinit = sol.u
    sol = solve_problem(m.second, f, D, rl2; u0=xinit)
end

function solve_problem(optimizer::BenchmarkSetup, f, D::Int, run_length::Int; u0 = pinit(D))

    method = optimizer.method

    optf = OptimizationFunction((u,_)->f(u), AutoForwardDiff())
    if optimizer.isboxed 
        prob = OptimizationProblem(optf, u0, lb = fill(-5.5, D), ub = fill(5.5, D))
    else
        prob = OptimizationProblem(optf, u0)
    end
    sol = solve(prob, method; maxiters = run_length)
    sol
end

function benchmark(
    optimizer::Union{Chain,BenchmarkSetup}, f::BBOBFunction, run_length::AbstractVector{Int}; 
    Ntrials::Int = 20, dimension::Int = 3, Δf::Real = 1e-6, CI_quantile=0.25, verbose=true
    )

    verbose && @info("$(string(optimizer))\t $f")

    t = (T) -> zeros(T, Ntrials, length(run_length))
        
    reached_minium = t(Bool)
    distance_to_xopt = t(Float64)
    fmin = t(Float64)
    callcount = t(Int)
    
    t = 0.0
    for j in 1:length(run_length)
        for i in 1:Ntrials
            try
                fcountner = FunctionCallsCounter(f.f)
                t += @elapsed sol = solve_problem(optimizer, fcountner, dimension, run_length[j])
                
                sol.objective, sol.u
                reached_minium[i,j] = sol.objective < Δf + f.f_opt
                fmin[i,j] = sol.objective - f.f_opt 
                distance_to_xopt[i,j] =  √sum(abs2.(sol.u - f.x_opt[1:dimension]))
                callcount[i,j] = fcountner.count

            catch err
                
                #throw(err)
                reached_minium[i,j] = false
                fmin[i,j] = NaN
                distance_to_xopt[i,j] = NaN
                callcount[i,j] = 0
                @warn(string(optimizer, " failed :", err))
            end
        end
    end
    t /= Ntrials*length(run_length)
    
    dr = x->dropdims(x, dims=1)
    success_rate = sum(reached_minium, dims=1)/Ntrials |> dr

    success_rate_qlow, success_rate_qhigh = compute_CI(success_rate, Ntrials, CI_quantile)
     
    BenchmarkResults(
        run_length = run_length,
        success_count = sum(reached_minium, dims=1) |> dr, 
        success_rate = success_rate, 
        success_rate_qlow = success_rate_qlow,
        success_rate_qhigh = success_rate_qhigh,
        distance_to_minimizer = mean(distance_to_xopt, dims=1) |> dr, 
        minimum = mean(fmin, dims=1) |> dr, 
        runtime = t,
        Neffective = Ntrials,
        callcount = mean(callcount, dims=1) |> dr,
        success_rate_per_function = [success_rate[end]]
    )
end


benchmark(optimizer, funcs, run_length::AbstractVector{Int}; 
    Ntrials::Int = 20, dimension::Int = 3, Δf::Real = 1e-6, CI_quantile=0.25
) = benchmark(
    BenchmarkSetup(optimizer), funcs, run_length; Ntrials, dimension, Δf, CI_quantile
)

#
function benchmark(
    optimizer::Union{Chain,BenchmarkSetup}, funcs::Vector{BBOBFunction}, run_length::AbstractVector{Int}; 
    Ntrials::Int = 20, dimension::Int = 3, Δf::Real = 1e-6, CI_quantile=0.25
    )
    
    res = [benchmark(optimizer, f, run_length; Ntrials, dimension, Δf) for f in funcs]
    reduce_res(res, field, f=mean) = f(getfield(r, field) for r in res)
    
    Neff = Ntrials * length(funcs)
    success_rate = reduce_res(res, :success_count, sum) / Neff
    
    success_rate_qlow, success_rate_qhigh = compute_CI(success_rate, Neff, CI_quantile)

    BenchmarkResults(
        run_length = run_length,
        success_count = reduce_res(res, :success_count, sum),
        success_rate = success_rate,
        success_rate_qlow = success_rate_qlow,
        success_rate_qhigh = success_rate_qhigh,
        distance_to_minimizer = reduce_res(res, :distance_to_minimizer),
        minimum = reduce_res(res, :minimum),
        runtime = reduce_res(res, :runtime, sum),
        Neffective = Neff,
        callcount = reduce_res(res, :callcount),
        success_rate_per_function = [res.success_rate[end] for res in res]
    )
end

##



