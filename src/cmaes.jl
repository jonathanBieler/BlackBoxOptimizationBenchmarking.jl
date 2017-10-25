
module CMAES

##

using Distributions, Compat, PDMats

@compat weights(μ) = normalize!([ (log.(μ+1) - log(i)) / (μ*log.(μ+1) -sum(log.(1:μ)) ) for i =1:μ],1)

function average_x!(out::Vector{T},x::Matrix{T},w::Vector{T},μ) where {T <:Number} 
    fill!(out,0)
    for i=1:μ
        for j=1:size(x,1)
            out[j] += x[j,i]*w[i] 
        end
    end
end

function init_parameters(xinit,λ,w,μ)
    n = length(xinit)
#    λ = round(Int, 3 + floor(3log(n)))
    
    μ_eff = 1.0 / sum(w.^2)
    c_σ =  (μ_eff + 2.0) / (n + μ_eff + 3.0)
    
    d_σ = 1.0 + 2.0*max(0, √((μ_eff-1)/(n+1)) -1) + c_σ
    
    c_c = 4.0/(n+4.0) 
    α = 2.0 / (n + √2)^2
    
    n, λ, μ_eff, c_σ, d_σ, c_c, α
end

function update_A_inv(A_inv, α, ν)

    1/√(1-α)*A_inv -1/√(1-α)/norm(ν)^2 *(
        1 - 1 / ( √(1 + α/(α/(1-α))*norm(ν)^2) )
    ) * ν * (ν' * A_inv)
end

function update_A(A, α, ν, p_c)

    √(1-α)*A + √(1-α)/norm(ν)^2 *(
        √(1 + α/(1-α)*norm(ν)^2) - 1
    ) * p_c * ν'
end

"""

Notation: α = c_cov

"""

function cmaes(f::Function, xinit::Vector{T}, σ::T, niter::Int, λ::Int) where T

    #init everything
    μ = floor(Int,λ/2)
    w = weights(μ)
    n, λ, μ_eff, c_σ, d_σ, c_c, α = init_parameters(xinit,λ,w,μ)
    
    A, A_inv = σ*eye(n), 1/σ*eye(n)
    p_σ, p_c = zeros(n), zeros(n)
    D = MultivariateNormal(zeros(n),PDMats.PDiagMat(ones(n)))
    
    #start looping
    
    x, z = zeros(n,λ), zeros(n,λ)
    x_w, z_w = zeros(n), zeros(n)
    
    c = zeros(λ,n)
    fx = zeros(λ)
    
    f_min, f_minp = Inf, Inf
    x_min, x_minp = zeros(n), zeros(n)
    tmp = zeros(n)

    for i = 1:niter
    
        # generate new points and evaluate fitness
        average_x!(x_w,x,w,μ)
        @inbounds for k = 1:λ
            z[:,k] .= rand(D)
            #x[:,k] .= x_w .+ σ*A*z[:,k]
            
            Base.BLAS.gemv!('N', σ, A, z[:,k], 0.0, tmp)#this isn't really faster
            for i=1:n
               x[i,k] = tmp[i] + x_w[i]
            end

            fx[k] = f(x[:,k])
        end
        
        # sort by fitness
        idx = sortperm(fx)
        x, z = x[:,idx], z[:,idx]
        
        # save best candidate
        if fx[idx[1]] < f_min
            f_minp = f_min
            x_minp = x_min

            x_min = x[:,1]
            f_min = fx[idx[1]]
        end

        # termination
        if i > 1

            f_tol = 1e-12
            x_tol = 1e-12
        
            abs(f_min - f_minp) / (abs(f_min) + f_tol) < f_tol && return (x_min,f_min)
            maximum( abs(x_min[i] - x_minp[i])  for i=1:n ) < x_tol && return (x_min,f_min)
        end
        
        # update
        average_x!(z_w,z,w,μ)
        
        p_c = (1 - c_c)*p_c + √(c_σ*(2 - c_σ)*μ_eff) * A * z_w
        
        ν = A_inv * p_c
        
        A_inv = update_A_inv(A_inv, α, ν)
        A = update_A(A, α, ν, p_c)
            
        p_σ = (1-c_σ) * p_σ  + √(c_σ*(2-c_σ)*μ_eff) * z_w
        
        chi_n = √n*(1-1/(4*n) + 1/(21*n^2))
        
        σ = σ * exp(c_σ/d_σ*(norm(p_σ)/chi_n -1))
    
    end

    x_min,f_min
end

##
end


## test
if false

using BBOBFunctions

D = 2
niter = 2000

x = randn(D)*10

f = BBOBFunctions.f3

@time pmin = CMAES.cmaes(f,x,1,niter,40)

f(BBOBFunctions.x3_opt) - BBOBFunctions.f3_opt

sleep(0.01)
println(f(pmin) - BBOBFunctions.f3_opt)
println(pmin - BBOBFunctions.x3_opt[1:D])

## compare to BlackBoxOptim


using BlackBoxOptim

fit = bboptimize(f, Method=:generating_set_search, SearchRange = (-15.0, 15.0), NumDimensions = D,MaxFuncEvals=niter)
pmin = best_candidate(fit)

sleep(0.1)
println(f(pmin) - BBOBFunctions.f3_opt)
println(pmin - BBOBFunctions.x3_opt[1:D])

## compare to other CMAES

include("cmaes2.jl")

pmin = cmaes(f,x,5*ones(D),stopeval=niter)

println(f(pmin) - BBOBFunctions.f3_opt)
println(pmin - BBOBFunctions.x3_opt[1:D])

## benchmark

Nrepeats = 30
niters = [100 500 1000 5000]
dimensions = [2]#collect(2:5:15)
Nmethod = 2
proportion = zeros(Nmethod,length(dimensions)*length(niters))
x_val = zeros(length(dimensions)*length(niters))

targetvalue = BBOBFunctions.f3_opt + 1/100*abs(BBOBFunctions.f3_opt)
f = BBOBFunctions.f3

ind_iter = 1
for i = 1:length(dimensions)
    
    D = dimensions[i]

    for niter in niters

        fvalues = zeros(Nmethod,Nrepeats)
        for r in 1:Nrepeats

            
            x = randn(D)*2
            pmin = CMAES.cmaes(f,x,1,niter,50)
            fvalues[1,r] = f(pmin)
            
            pmin = best_candidate(
                bboptimize(f, Method=:generating_set_search, SearchRange = (-15.0, 15.0), NumDimensions = D,MaxFuncEvals=niter)
            )
            fvalues[2,r] = f(pmin)

        end

        warn(ind_iter)

        for k=1:size(proportion,1)
            proportion[k,ind_iter] = mean(fvalues[k,:] .< targetvalue)
        end
        x_val[ind_iter] = niter/D
        ind_iter += 1
        
    end
end

p = plot(
    layer(x=x_val,y=proportion[1,:],Geom.point,Geom.line),
    layer(x=x_val,y=proportion[2,:],Geom.point,Geom.line,col("red"))
)

end
##





