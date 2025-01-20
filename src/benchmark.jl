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
end
BenchmarkResults(;run_length, success_count, success_rate, success_rate_qlow, success_rate_qhigh, distance_to_minimizer, minimum, runtime, Neffective) =
    BenchmarkResults(run_length, success_count, success_rate, success_rate_qlow, success_rate_qhigh, distance_to_minimizer, minimum, runtime, Neffective)


show(io::IO, b::BenchmarkResults) =  begin
    println(io, "BenchmarkResults :")
    print(io, "Run length : ")
    show(IOContext(io, :limit => true, :compact => true), b.run_length)
    print(io, "\nSuccess rate : ")
    show(IOContext(io, :limit => true, :compact => true), b.success_rate)
end

# https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Clopper%E2%80%93Pearson_interval
compute_CI(rate::Real, N::Int, q::Real) = quantile(Beta(1 + round(Int, N*rate), 1 + N - round(Int, N*rate)), q)

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

function benchmark(optimizer, f::BBOBFunction, run_length::AbstractVector{Int}; Ntrials::Int = 20, dimension::Int = 3, Δf::Real = 1e-6, CI_quantile=0.25)

    @info("$(string(optimizer))\t $f")

    t = (T) -> zeros(T, Ntrials, length(run_length))
        
    reached_minium = t(Bool)
    distance_to_xopt = t(Float64)
    fmin = t(Float64)
    
    t = 0.0
    for j in 1:length(run_length)
        for i in 1:Ntrials
            try
                t += @elapsed mfit = optimize(optimizer, f.f, dimension, run_length[j])
                
                reached_minium[i, j] = minimum(mfit) < Δf + f.f_opt
                fmin[i, j] = minimum(mfit) - f.f_opt 
                distance_to_xopt[i, j] =  √sum(abs2.(minimizer(mfit) - f.x_opt[1:dimension]))
                
            catch err
                @warn(err)
                
                reached_minium[i, j] = false
                fmin[i, j] = NaN
                distance_to_xopt[i, j] = NaN
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
        Neffective = Ntrials
    )
end

#
function benchmark(optimizer, funcs::Vector{BBOBFunction}, run_length::AbstractVector{Int}; Ntrials::Int = 20, dimension::Int = 3, Δf::Real = 1e-6, CI_quantile=0.25)
    
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
        Neffective = Neff
    )
end




##



