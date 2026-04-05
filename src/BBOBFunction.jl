## Constants and Functions

using StaticArrays, LinearAlgebra, Random

const maximum_dimension = 100

@inline _safe_sqrt(x) = sqrt(max(x, zero(x)))

@inline function T_osz(xi::T) where T <: Number
    xhat = xi != zero(T) ? log(abs(xi)) : zero(T)
    c1 = xi > zero(T) ? T(10) : T(5.5)
    c2 = xi > zero(T) ? T(7.9) : T(3.1)
    sign(xi) * exp(xhat + T(0.049) * (sin(c1 * xhat) + sin(c2 * xhat)))
end

@inline T_osz(x::SVector) = map(T_osz, x)

@inline function T_asy(x::SVector{N, T}, β) where {N, T}
    SVector{N, T}(ntuple(Val(N)) do i
        xi = x[i]
        xi > zero(T) ? xi^(one(T) + T(β) * T(i - 1) / T(N - 1) * _safe_sqrt(xi)) : xi
    end)
end

@inline function Λ_mul(::Val{N}, α::T, x::SVector{N, T}) where {N, T}
    SVector{N, T}(ntuple(i -> α^(T(0.5) * T(i - 1) / T(N - 1)) * x[i], Val(N)))
end

@inline function f_pen(x::SVector{N, T}) where {N, T}
    sum(max.(zero(T), abs.(x) .- T(5)) .^ 2)
end

@inline function ellip_weights(::Val{N}, ::Type{T}, exponent::T) where {N, T}
    SVector{N, T}(ntuple(i -> T(10)^(exponent * T(i - 1) / T(N - 1)), Val(N)))
end

## BBOBFunction

struct BBOBFunction{F, N, M}
    name::String
    f::F
    x_opt::SVector{N, Float32}
    f_opt::Float32
    Q::SMatrix{N, N, Float32, M}
    R::SMatrix{N, N, Float32, M}
end

function (func::BBOBFunction{F, N, M})(x) where {F, N, M}
    x_static = SVector{N, Float32}(x)
    Float64(func.f(x_static, func.x_opt, func.f_opt, func.Q, func.R))
end

Base.show(io::IO, f::BBOBFunction) = print(io, f.name)
Base.broadcastable(f::BBOBFunction) = Ref(f)
minimum(f::BBOBFunction) = f.f_opt
minimizer(f::BBOBFunction) = f.x_opt

function make_rotation(::Val{N}, seed::Int) where N
    rng = MersenneTwister(seed)
    A = randn(rng, Float32, N, N)
    SMatrix{N, N, Float32}(Matrix(qr(A).Q))
end

function make_x_opt(::Val{N}, seed::Int) where N
    rng = MersenneTwister(seed)
    SVector{N, Float32}(rand(rng, Float32, N) .* 10f0 .- 5f0)
end

function make_x_opt_linear_slope(::Val{N}, seed::Int) where N
    rng = MersenneTwister(seed)
    SVector{N, Float32}(ntuple(i -> rand(rng) < 0.5f0 ? -5f0 : 5f0, Val(N)))
end

function make_x_opt_schwefel(::Val{N}, seed::Int) where N
    rng = MersenneTwister(seed)
    val = Float32(4.2096874633 / 2)
    SVector{N, Float32}(ntuple(i -> rand(rng) < 0.5f0 ? -val : val, Val(N)))
end

function make_f_opt(seed::Int)
    rng = MersenneTwister(seed)
    Float32(clamp(round(randn(rng) * 100, digits = 2), -1000, 1000))
end

## f1, Sphere Function

""" Sphere Function """
@inline function f1(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = x .- x_opt
    sum(z .^ 2) + f_opt
end

## f2, Ellipsoidal Function

""" Ellipsoidal Function """
@inline function f2(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = T_osz(x .- x_opt)
    w = ellip_weights(Val(N), T, T(6))
    sum(w .* z .^ 2) + f_opt
end

## f3, Rastrigin Function

""" Rastrigin Function """
@inline function f3(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = Λ_mul(Val(N), T(10), T_asy(T_osz(x .- x_opt), T(0.2)))
    T(10) * (T(N) - sum(cos.(T(2π) .* z))) + sum(z .^ 2) + f_opt
end

## f4, Buche-Rastrigin Function

""" Buche-Rastrigin Function """
@inline function f4(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = T_osz(x .- x_opt)
    s = SVector{N, T}(ntuple(i -> isodd(i) ? T(10) * T(10)^(T(0.5) * T(i - 1) / T(N - 1)) :
                                              T(10)^(T(0.5) * T(i - 1) / T(N - 1)), Val(N)))
    z = s .* z
    T(10) * (T(N) - sum(cos.(T(2π) .* z))) + sum(z .^ 2) + T(100) * f_pen(x) + f_opt
end

## f5, Linear Slope

""" Linear Slope """
@inline function f5(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    s = SVector{N, T}(ntuple(i -> sign(x_opt[i]) * T(10)^(T(i - 1) / T(N - 1)), Val(N)))
    z = SVector{N, T}(ntuple(i -> x_opt[i] * x[i] < T(25) ? x[i] : x_opt[i], Val(N)))
    sum(T(5) .* abs.(s) .- s .* z) + f_opt
end

## f6, Attractive Sector Function

""" Attractive Sector Function """
@inline function f6(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = Q * Λ_mul(Val(N), T(10), R * (x .- x_opt))
    z = SVector{N, T}(ntuple(i -> x_opt[i] * z[i] > zero(T) ? T(100) * z[i] : z[i], Val(N)))
    T_osz(sum(z .^ 2))^T(0.9) + f_opt
end

## f7, Step Ellipsoidal Function

""" Step Ellipsoidal Function """
@inline function f7(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = Λ_mul(Val(N), T(10), R * (x .- x_opt))
    zhat_1 = z[1]
    z = SVector{N, T}(ntuple(i -> z[i] > T(0.5) ?
        floor(T(0.5) + z[i]) :
        floor(T(0.5) + T(10) * z[i]) / T(10), Val(N)))
    z = Q * z
    w = ellip_weights(Val(N), T, T(2))
    T(0.1) * max(abs(zhat_1) / T(1e4), sum(w .* z .^ 2)) + f_pen(x) + f_opt
end

## f8, Rosenbrock Function, original

""" Rosenbrock Function, original """
@inline function f8(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = max(one(T), T(sqrt(N)) / T(8)) .* (x .- x_opt) .+ one(T)
    v = SVector{N - 1, T}(ntuple(i -> T(100) * (z[i]^2 - z[i + 1])^2 + (z[i] - one(T))^2, Val(N - 1)))
    sum(v) + f_opt
end

## f9, Rosenbrock Function, rotated

""" Rosenbrock Function, rotated """
@inline function f9(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = max(one(T), T(sqrt(N)) / T(8)) .* (R * (x .- x_opt)) .+ one(T)
    v = SVector{N - 1, T}(ntuple(i -> T(100) * (z[i]^2 - z[i + 1])^2 + (z[i] - one(T))^2, Val(N - 1)))
    sum(v) + f_opt
end

## f10, Ellipsoidal Function

""" Ellipsoidal Function """
@inline function f10(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = T_osz(R * (x .- x_opt))
    w = ellip_weights(Val(N), T, T(2))
    sum(w .* z .^ 2) + f_opt
end

## f11, Discus Function

""" Discus Function """
@inline function f11(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = T_osz(R * (x .- x_opt))
    T(1e6) * z[1]^2 + sum(z .^ 2) - z[1]^2 + f_opt
end

## f12, Bent Cigar Function

""" Bent Cigar Function """
@inline function f12(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = R * T_asy(R * (x .- x_opt), T(0.5))
    z[1]^2 + T(1e6) * (sum(z .^ 2) - z[1]^2) + f_opt
end

## f13, Sharp Ridge Function

""" Sharp Ridge Function """
@inline function f13(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = Q * Λ_mul(Val(N), T(10), R * (x .- x_opt))
    z[1]^2 + T(100) * _safe_sqrt(sum(z .^ 2) - z[1]^2) + f_opt
end

## f14, Different Powers Function

""" Different Powers Function """
@inline function f14(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = R * (x .- x_opt)
    pw = SVector{N, T}(ntuple(i -> abs(z[i])^(T(2) + T(4) * T(i - 1) / T(N - 1)), Val(N)))
    _safe_sqrt(sum(pw)) + f_opt
end

## f15, Rastrigin Function

""" Rastrigin Function """
@inline function f15(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = R * Λ_mul(Val(N), T(10), Q * T_asy(T_osz(R * (x .- x_opt)), T(0.2)))
    T(10) * (T(N) - sum(cos.(T(2π) .* z))) + sum(z .^ 2) + f_opt
end

## f16, Weierstrass Function

""" Weierstrass Function """
@inline function f16(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = R * Λ_mul(Val(N), T(1 / 100), Q * T_osz(R * (x .- x_opt)))
    f0 = zero(T)
    for k in 0:11
        f0 += T(1) / T(2)^k * cos(T(2π) * T(3)^k * T(0.5))
    end
    s = zero(T)
    for j in 1:N
        for k in 0:11
            s += T(1) / T(2)^k * cos(T(2π) * T(3)^k * (z[j] + T(0.5)))
        end
    end
    T(10) * (T(1) / T(N) * s - f0)^3 + T(10) / T(N) * f_pen(x) + f_opt
end

## f17, Schaffers F7 Function

""" Schaffers F7 Function """
@inline function f17(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = Λ_mul(Val(N), T(10), Q * T_asy(R * (x .- x_opt), T(0.5)))
    s = SVector{N - 1, T}(ntuple(i -> _safe_sqrt(z[i]^2 + z[i + 1]^2), Val(N - 1)))
    v = SVector{N - 1, T}(ntuple(i -> _safe_sqrt(s[i]) * (T(1) + sin(T(50) * s[i]^T(0.2))^2), Val(N - 1)))
    (T(1) / T(N - 1) * sum(v))^2 + T(10) * f_pen(x) + f_opt
end

## f18, Schaffers F7 Function, moderately ill-conditioned

""" Schaffers F7 Function, moderately ill-conditioned """
@inline function f18(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = Λ_mul(Val(N), T(1000), Q * T_asy(R * (x .- x_opt), T(0.5)))
    s = SVector{N - 1, T}(ntuple(i -> _safe_sqrt(z[i]^2 + z[i + 1]^2), Val(N - 1)))
    v = SVector{N - 1, T}(ntuple(i -> _safe_sqrt(s[i]) * (T(1) + sin(T(50) * s[i]^T(0.2))^2), Val(N - 1)))
    (T(1) / T(N - 1) * sum(v))^2 + T(10) * f_pen(x) + f_opt
end

## f19, Composite Griewank-Rosenbrock Function F8F2

""" Composite Griewank-Rosenbrock Function F8F2 """
@inline function f19(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    z = max(one(T), T(sqrt(N)) / T(8)) .* (R * (x .- x_opt)) .+ one(T)
    s = SVector{N - 1, T}(ntuple(i -> T(100) * (z[i]^2 - z[i + 1])^2 + (z[i] - one(T))^2, Val(N - 1)))
    v = SVector{N - 1, T}(ntuple(i -> s[i] / T(4000) - cos(s[i]), Val(N - 1)))
    T(10) / T(N - 1) * sum(v) + T(10) + f_opt
end

## f20, Schwefel Function

""" Schwefel Function """
@inline function f20(x::SVector{N, T}, x_opt, f_opt, Q, R) where {N, T}
    one_pm = sign.(x_opt)
    x_scaled = T(2) .* one_pm .* x

    abs_x_opt = abs.(x_opt)

    z = SVector{N, T}(ntuple(Val(N)) do i
        i == 1 ? x_scaled[1] : x_scaled[i] + T(0.25) * (x_scaled[i - 1] - T(2) * abs_x_opt[i - 1])
    end)

    z_shifted = z .- T(2) .* abs_x_opt
    z = T(100) .* (Λ_mul(Val(N), T(10), z_shifted) .+ T(2) .* abs_x_opt)

    v = SVector{N, T}(ntuple(i -> z[i] * sin(_safe_sqrt(abs(z[i]))), Val(N)))
    -T(1) / (T(100) * T(N)) * sum(v) + T(4.189828872724339) + T(100) * f_pen(z ./ T(100)) + f_opt
end

const BBOB_FUNCTIONS = [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10,
                        f11, f12, f13, f14, f15, f16, f17, f18, f19, f20]

const BBOB_NAMES = [
    "F1  Sphere", "F2  Ellipsoidal", "F3  Rastrigin", "F4  Buche-Rastrigin",
    "F5  Linear Slope", "F6  Attractive Sector", "F7  Step Ellipsoidal",
    "F8  Rosenbrock", "F9  Rosenbrock Rotated", "F10 Ellipsoidal 2",
    "F11 Discus", "F12 Bent Cigar", "F13 Sharp Ridge", "F14 Different Powers",
    "F15 Rastrigin 2", "F16 Weierstrass", "F17 Schaffers F7",
    "F18 Schaffers F7 Ill-Cond", "F19 Griewank-Rosenbrock", "F20 Schwefel"]

"""
    bbob_suite(Val(N); seed=42) -> Vector{BBOBFunction}

Build the full 20-function BBOB suite for dimension `N`.
Each function gets deterministic random rotations and optima from `seed`.
"""
function bbob_suite(::Val{N}; seed = 42) where N
    suite = BBOBFunction[]
    for (i, fn) in enumerate(BBOB_FUNCTIONS)
        if fn === f5
            x_opt = make_x_opt_linear_slope(Val(N), seed + i)
        elseif fn === f20
            x_opt = make_x_opt_schwefel(Val(N), seed + i)
        else
            x_opt = make_x_opt(Val(N), seed + i)
        end
        f_opt = make_f_opt(seed + 100 + i)
        Qmat = make_rotation(Val(N), seed + 200 + i)
        Rmat = make_rotation(Val(N), seed + 300 + i)
        push!(suite, BBOBFunction(BBOB_NAMES[i], fn, x_opt, f_opt, Qmat, Rmat))
    end
    suite
end

list_functions() = BBOB_NAMES