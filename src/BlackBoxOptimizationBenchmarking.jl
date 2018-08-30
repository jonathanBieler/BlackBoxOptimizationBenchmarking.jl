module BlackBoxOptimizationBenchmarking

    using Distributions, Memoize, Compat, Optim
    using LinearAlgebra, Distributed

    import Base: enumerate, show, string
    export show, BBOBFunction, enumerate
    
## Constants and Functions

    const maximum_dimension = 100

    T_asy(x::Vector{T},β) where T <: Number = [x[i] > 0 ? x[i]^(1+β*(i-1)/(length(x)-1)*√x[i]) : x[i] for i=1:length(x)]

    function T_asy(x::Matrix{T},β) where T <: Number 
        z = similar(x)
        for i=1:size(x,1)
            z[i,:] = T_asy(x[i,:],β)
        end
        z
    end
    
    function T_osz(xi::T) where T <: Number 
        xhat = xi != zero(T) ? log(abs(xi)) : zero(T)
        c1 = xi > 0 ? 10. : 5.5
        c2 = xi > 0 ? 7.9 : 3.1

        sign(xi) * exp(xhat + 0.049*(sin(c1*xhat)+sin(c2*xhat) ))
    end

    function T_osz(x::Vector{T}) where T <: Number 
        z = similar(x)
        for i=1:length(x)
            z[i] = T_osz(x[i])
        end
        z
    end

    function T_osz(x::Matrix{T}) where T <: Number 
        z = similar(x)
        for i=1:size(x,1)
            z[i,:] = T_osz(x[i,:])
        end
        z
    end

    function Λ(α,D)
        m = zeros(D,D)
        for i=1:D
            m[i,i] = α^(1/2 *(i-1)/(D-1))
        end
        m
    end
    
    C(i,D,n) = 10^(2*(i-1)/(D-1))#conditioning
    
    f_pen(x) =  sum(max(0,abs(xi)-5)^2 for xi in x)
    
    const one_pm = sign.(randn(maximum_dimension))
    
    #rotation matrices (probably a bit wrong)
    @memoize function Q(D)
        r = randn(D); r = r/norm(r)
        Q = [r nullspace(Matrix(r'))]
    end
    @memoize function R(D)
        r = randn(D); r = r/norm(r)
        R = [r nullspace(Matrix(r'))]
    end
    
    ∑(x) = sum(x)

## BBOBFunction

    struct BBOBFunction{F<:Function}
        name::String
        f::F
        x_opt::Array{Float64,1}
        f_opt::Float64
    end

    (f::BBOBFunction)(x) = f.f(x)
    show(io::IO,f::BBOBFunction) =  print(io,f.name)
    Base.broadcastable(f::BBOBFunction) = Ref(f)

    test_x_opt(f::BBOBFunction) = @assert f(f.x_opt) ≈ f.f_opt
    
## helpers to define function

    fun_symbols(n) = (map(Symbol, ["F$(n)", "f$(n)", "x$(n)_opt", "f$(n)_opt"])...,)

    macro BBOBFunction(name,n)
        F,f,x_opt,f_opt = fun_symbols(n)
        esc(:( $F = BBOBFunction($name,$f,$x_opt,$f_opt) ))
    end

    function enumerate(::Type{BBOBFunction})
         n = names(BlackBoxOptimizationBenchmarking, all=true)
         v = map(x->getfield(BlackBoxOptimizationBenchmarking,x),n)
         
         out = BBOBFunction[]
         for x in v 
             typeof(x) <: BBOBFunction && push!(out,x)
         end
         out
    end

    """
        Define x1_opt and f1_opt.

    """
    macro define_x_and_f_opt(n)
        F,f,x_opt,f_opt = fun_symbols(n)
        esc(quote
            const $x_opt = rand(Uniform(-4,4),maximum_dimension)
            const $f_opt = min(1000,max(-1000, round(rand(Cauchy(0,100)), digits=2)))
        end)
    end
    
## Functions

    ## f1, Sphere Function

    @define_x_and_f_opt(1)

    """ Sphere Function """
    f1(x) = ∑( (x[i]-x1_opt[i])^2 for i=1:length(x) ) + f1_opt

    @BBOBFunction("Sphere",1)

    ## f2, Ellipsoidal Function

    @define_x_and_f_opt(2)

    """ Ellipsoidal Function """
    function f2(x) 
        D = length(x)
        z = T_osz(x-x2_opt[1:D])
        ∑( 10^(6*(i-1)/(D-1)) * z[i]^2 for i=1:length(x) ) + f2_opt
    end

    @BBOBFunction("Ellipsoidal",2)

    ## f3, Rastrigin Function

    @define_x_and_f_opt(3)

    """ Rastrigin Function """
    function f3(x) 
        D = length(x)
        z = Λ(10,D) * T_asy(T_osz(x-x3_opt[1:D]),0.2)
        10*(D - ∑( cos(2π*z[i]) for i=1:D )) + norm(z)^2 + f3_opt
    end

    @BBOBFunction("Rastrigin",3)
    
    ## f4, Buche-Rastrigin Function
    
    @define_x_and_f_opt(4)

    """ Buche-Rastrigin Function """
    function f4(x) 
        D = length(x)
        z = T_osz(x-x4_opt[1:D])
        s = [isodd(i) ? 10*10^(0.5*(i-1)/(D-1)) : 10^(0.5*(i-1)/(D-1)) for i=1:D] 

        for i=1:D 
            @inbounds z[i] = s[i]*z[i] 
        end

        10*(D - ∑( cos(2π*z[i]) for i=1:D )) + norm(z)^2 + 100*f_pen(x) + f4_opt
    end

    @BBOBFunction("Buche-Rastrigin",4)
    
    ## f5, Linear Slope
    
    const x5_opt = 5*one_pm
    const f5_opt = min(1000,max(-1000, round(rand(Cauchy(0,100)), digits=2)))

    """ Linear Slope """
    function f5(x) 
        D = length(x)
        z = [ x5_opt[i]*x[i] < 25 ? x[i] : x5_opt[i] for i=1:D ]
        s = [sign(x5_opt[i])*10^((i-1)/(D-1))  for i=1:D] 

        ∑( 5*abs(s[i]) -s[i]*z[i] for i=1:D ) + f5_opt
    end

    @BBOBFunction("Linear Slope",5)
    
    ## f6, Attractive Sector Function
    
    @define_x_and_f_opt(6)

    """ Attractive Sector Function """
    function f6(x) 
        D = length(x)
        
        z = Q(D)*Λ(10,D)*R(D)*(x - x6_opt[1:D])
        
        @inbounds for i=1:D 
            z[i] = x6_opt[i]*z[i] > 0 ? 100*z[i] : z[i]
        end
        
        T_osz( ∑( z[i]^2 for i=1:D ))^0.9 + f6_opt
    end

    @BBOBFunction("Attractive Sector",6)

    
    ## f7, Step Ellipsoidal Function
    
    @define_x_and_f_opt(7)

    """ Step Ellipsoidal Function """
    function f7(x) 
        D = length(x)
        
        z = Λ(10,D)*R(D)*(x - x7_opt[1:D])
        zhat_1 = copy(z[1])
        
        @inbounds for i=1:D 
            z[i] = z[i] > 0.5 ? floor(0.5 + z[i]) : floor(0.5 + 10*z[i])/10
        end
        z = Q(D)*z
         
        0.1*max(abs(zhat_1)/(10^4), ∑( 10^(2*(i-1)/(D-1)) * z[i]^2 for i=1:D ) ) + f_pen(x) + f7_opt
    end

    @BBOBFunction("Step Ellipsoidal Function",7)
    
    ## f8, Rosenbrock Function, original
    
    @define_x_and_f_opt(8)

    """ Rosenbrock Function, original """
    function f8(x) 
        D = length(x)
        
        z = max(1,√D/8)*(x - x8_opt[1:D]) .+ 1
        
        ∑( 100*(z[i]^2 - z[i+1])^2 + (z[i]-1)^2 for i=1:D-1 ) + f8_opt
    end

    @BBOBFunction("Rosenbrock Function, original",8)
    
    ## f9, Rosenbrock Function, rotated
    
    @define_x_and_f_opt(9)

    """ Rosenbrock Function, rotated """
    function f9(x) 
        D = length(x)
        
        z = max(1,√D/8)*R(D)*(x - x9_opt[1:D]) .+ 1
        
        ∑( 100*(z[i]^2 - z[i+1])^2 + (z[i]-1)^2 for i=1:D-1 ) + f9_opt
    end

    @BBOBFunction("Rosenbrock Function, rotated",9)
    
    ## f10, Ellipsoidal Function
    
    @define_x_and_f_opt(10)

    """ Ellipsoidal Function """
    function f10(x) 
        D = length(x)
        
        z = T_osz( R(D)*(x-x10_opt[1:D])) 
        
        ∑( C(i,D,6)*z[i]^2 for i=1:D ) + f10_opt
    end

    @BBOBFunction("Ellipsoidal Function",10)
    
    ## f11, Discus Function
    
    @define_x_and_f_opt(11)

    """ Discus Function """
    function f11(x) 
        D = length(x)
        
        z = T_osz( R(D)*(x - x11_opt[1:D])) 
        
        10^6*z[1]^2 + ∑( z[i]^2 for i=2:D ) + f11_opt
    end

    @BBOBFunction("Discus Function",11)
    
    
    ## f12, Bent Cigar Function
    
    @define_x_and_f_opt(12)

    """ Bent Cigar Function """
    function f12(x) 
        D = length(x)
        
        z = R(D)*T_asy( R(D)*(x - x12_opt[1:D]), 0.5) 
        
        z[1]^2 + 10^6*∑( z[i]^2 for i=2:D ) + f12_opt
    end

    @BBOBFunction("Bent Cigar Function",12)
    
    
    ## f13, Sharp Ridge Function
    
    @define_x_and_f_opt(13)

    """ Sharp Ridge Function """
    function f13(x)
        D = length(x)
        
        z = Q(D)*Λ(10,D)*R(D)*(x - x13_opt[1:D])
        
        z[1]^2 + 100*√(∑(z[i]^2 for i=2:D )) + f13_opt
    end

    @BBOBFunction("Sharp Ridge Function",13)
    
    
    ## f14, Different Powers Function
    
    @define_x_and_f_opt(14)

    """ Different Powers Function """
    function f14(x)
        D = length(x)
        
        z = R(D)*(x - x14_opt[1:D])
        
        √(∑(abs(z[i])^(2+4(i-1)/(D-1)) for i=1:D )) + f14_opt
    end

    @BBOBFunction("Different Powers Function",14)
    
    ## f15, Rastrigin Function

    @define_x_and_f_opt(15)
    
    """ Rastrigin Function """
    function f15(x)
        D = length(x)

        z = R(D)*Λ(10,D)*Q(D)*T_asy(T_osz( R(D)*(x-x15_opt[1:D]) ),0.2)
        10*(D - ∑( cos(2π*z[i]) for i=1:D )) + norm(z)^2 + f15_opt
    end
    
    @BBOBFunction("Rastrigin Function",15)

    ## f16, Weierstrass Function

    @define_x_and_f_opt(16)
    
    const f0 = ∑( 1/2^k * cos(2π*3^k*1/2) for k=0:11 )

    """ Weierstrass Function """
    function f16(x)
        D = length(x)

        z = R(D)*Λ(1/100,D)*Q(D)*T_osz( R(D)*(x-x16_opt[1:D]) )
        s = ∑( ∑( 1/2^k * cos(2π*3^k*(z[i]+1/2)) for k=0:11 ) for i=1:D )
        10*( 1/D*s - f0 )^3 + 10/D*f_pen(x) + f16_opt
    end
    
    @BBOBFunction("Weierstrass Function",16)

    ## f17, Schaffers F7 Function
    
    @define_x_and_f_opt(17)
    
    """ Schaffers F7 Function"""
    function f17(x)
        D = length(x)

        z = Λ(10,D)*Q(D)*T_asy( R(D)*(x-x17_opt[1:D]), 0.5)
        s = [√(z[i]^2 + z[i+1]^2) for i=1:D-1]
        
        (1/(D-1)*∑( √s[i]*(1 + sin(50*s[i]^1/5 )^2 ) for i=1:D-1 ))^2 + 10*f_pen(x) + f17_opt
    end
    
    @BBOBFunction("Schaffers F7 Function",17)

    ## f18, Schaffers F7 Function, moderately ill-conditioned

    @define_x_and_f_opt(18)
    
    """ Schaffers F7 Function, moderately ill-conditioned """
    function f18(x)
        D = length(x)

        z = Λ(1000,D)*Q(D)*T_asy( R(D)*(x-x18_opt[1:D]), 0.5)
        s = [√(z[i]^2 + z[i+1]^2) for i=1:D-1]
        
        (1/(D-1)*∑( √s[i]*(1 + sin(50*s[i]^1/5 )^2 ) for i=1:D-1 ))^2 + 10*f_pen(x) + f18_opt
    end
    
    @BBOBFunction("Schaffers F7 Function, moderately ill-conditioned",18)
        

    ## f19, Composite Griewank-Rosenbrock Function F8F2

    @define_x_and_f_opt(19)
    
    """ Composite Griewank-Rosenbrock Function F8F2 """
    function f19(x)
        D = length(x)

        z = max(1,√D/8)*R(D)*(x - x19_opt[1:D]) .+ 1
        s = [ 100*(z[i]^2 - z[i+1])^2 + (z[i]-1)^2 for i=1:D-1]
        
        10/(D-1)*∑( s[i]/4000 -cos(s[i]) for i=1:D-1 ) + 10 + f19_opt
    end
    
    @BBOBFunction("Composite Griewank-Rosenbrock Function F8F2",19)


    ## f20, Schwefel Function
    # https://github.com/numbbo/coco/issues/837

    const x20_opt = 4.2096874633/2 * one_pm
    const f20_opt = min(1000,max(-1000, round(rand(Cauchy(0,100)), digits=2)))
    
    """ Schwefel Function """
    function f20(x)
        D = length(x)

        x = 2*one_pm[1:D] .* x
        
        z = [i==1 ? x[1] : x[i] + 0.25*(x[i-1]-2*abs(x20_opt[i-1]) ) for i=1:D]
        z = 100*( Λ(10,D)*(z- 2*abs.(x20_opt[1:D]) ) + 2*abs.(x20_opt[1:D]) )
        
        -1/(100*D)*∑( z[i]*sin(√(abs(z[i]))) for i=1:D ) + 4.189828872724339 + 100*f_pen(z/100) + f20_opt
    end
    
    @BBOBFunction("Schwefel Function",20)

## Tests

    map(test_x_opt,enumerate(BBOBFunction))

## utilities 
    
    include("benchmark.jl")

    function run_benchmark()
        include(joinpath(Pkg.dir(),"BlackBoxOptimizationBenchmarking","src","run_benchmark.jl"))
    end

    #x = [collect(linspace(-10,10,500)) collect(linspace(-10,10,500))]
    #
    #T_osz(x)
    #
    #y = [f3(x[i,:]) for i=1:size(x,1)]
    #
    #plot(x=x[:,1],y=y,Geom.line)

##
end # module
