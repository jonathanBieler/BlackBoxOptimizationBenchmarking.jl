##

function benchmark(optimizer, f::BBOBFunction, run_lengths, Ntrials, dimensions, Δf)

    @info("$(string(optimizer))\t $f")

    t = (T) -> NamedDimsArray{(:trial, :run_length, :dimension)}(
        zeros(T, Ntrials, length(run_lengths), length(dimensions))
    )

    reached_minium = t(Bool)
    distance_to_xopt = t(Float64)
    fmin = t(Float64)
    
    t = 0.0
    for j in 1:length(run_lengths), k in 1:length(dimensions)
        for i in 1:Ntrials
            try
                D = dimensions[k]
                t += @elapsed mfit = optimize(optimizer, f.f, D, run_lengths[j])
                
                reached_minium[i, j, k] = minimum(mfit) < Δf + f.f_opt
                fmin[i, j, k] = minimum(mfit) - f.f_opt 
                distance_to_xopt[i, j, k] =  √sum(abs2.(minimizer(mfit) - f.x_opt[1:D]))
                
            catch err
                @warn(err)
                
                reached_minium[i, j, k] = false
                fmin[i, j, k] = NaN
                distance_to_xopt[i, j, k] = NaN
                @warn(string(optimizer, " failed :", err))
            end
        end
    end
    t /= Ntrials*length(run_lengths)*length(dimensions)
    
    (
     success_counts = sum(reached_minium, dims=1), 
     distance_to_minimizer = mean(distance_to_xopt, dims=1), 
     minimum = mean(fmin, dims=1), 
     runtime = t
    )
end

#
function benchmark(optimizer, funcs::Vector, run_lengths, Ntrials, dimensions, Δf)
    
    res = [benchmark(optimizer, f, run_lengths, Ntrials, dimensions, Δf) for f in funcs]
    reduce_res(res, field, f=mean) = f(getfield(r, field) for r in res)
    
    (
        success_counts = reduce_res(res, :success_counts, sum),
        distance_to_minimizer = reduce_res(res, :distance_to_minimizer),
        minimum = reduce_res(res, :minimum),
        runtime = reduce_res(res, :runtime, sum),
    )
end




##



