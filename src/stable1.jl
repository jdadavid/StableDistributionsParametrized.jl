# extension to StableDistributions to using Nolan's parametrizations (0,1,2)
# knowing that StableDistributions implicit parametrization is 1, whereas Nolan uses 0 by default in his book

# methods mainly from John P. Nolan, "Univariate Stable Distributions", Springer 2020

"""
    Stable1(α, β, γ, δ)

The *Stable distribution* with stability index 0 < `α` ≤ 2, skewness parameter -1 ≤ `β` ≤ 1, scale 0 < `γ` and location `δ`
with Nolan'sparametrization 1

From from John P. Nolan, "Univariate Stable Distributions", Springer 2020, chapter 1.3 (p.7)

    p=0  is best for numerical works and statisticsinference; the characteristic function cf is continuous in all parameters

    p=1  is best for working with characteristic function and nice algebraic properties
         but location of the mode is unbounded near α=1 (tends to +Inf if α<1 and to -Inf if α>1)

    p=2  a.k.a mode-centered, is best for users
         location and scale parameters agree with standard Normal's when α=2, and with Cauchy's when α=1

```julia

Stable1(α)            # standard symmetric α-Stable distribution equivalent to Stable1(α, 0, 1, 0)
Stable1(α, β)         # standard α-Stable distribution with skewness parameter β equivalent to Stable1(α, β, 1, 0)
Stable1(α, β, γ, δ)   # α-Stable distribution with skewness parameter β, scale γ and location δ

params(d)             # Get the parameters, i.e. (α, β, γ, δ)
shape(d)              # Get the shape, i.e. (α, β)
location(d)           # Get the location, i.e. δ
scale(d)              # Get the scale, i.e. γ
minimum(d)            # Get the lower bound
maximum(d)            # Get the upper bound
mode(d)               # Get the mode (location of pdf maximum)
```

External links

* [Stable distribution on Wikipedia](https://en.wikipedia.org/wiki/Stable_distribution)

"""
struct Stable1{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    β::T
    γ::T
    δ::T
    α1::T
    β1::T
    σ1::T
    μ1::T
    function Stable1{T}(α, β, γ, δ; p, check_args::Bool=true) where {T}
        @check_args Stable (α, zero(T) < α <= 2one(T)) (β, -one(T) <= β <= one(T)) (γ, zero(T) < γ) ((α, β), α != 2 || β == 0 )
        (α1,β1,σ1,μ1) = _paramconvstabledists(2,1,α,β,γ,δ)
        new{T}(α, β, γ, δ, α1, β1, σ1, μ1)
    end
end


StablefromStable1(d1::Stable1)   = Stable(params(d1)...)
Stable1fromStable(ds::Stable)     = Stable1(params(ds)...)
StablepfromStable1(d1::Stable1)  = Stablep(params(d1)...; p=1)
Stable1fromStablep(dp::Stablep)  = dp.p == 1 ? Stable1(params(dp)...) : nothing

Stable1(α::T, β::T, γ::T, δ::T; check_args::Bool=true) where {T <: Real} = Stable1{T}(α, β, γ, δ; check_args=check_args)

Stable1(α::Real, β::Real, γ::Real, δ::Real; check_args::Bool=true) = Stable1(promote(α, β, γ, δ)...; check_args=check_args)
Stable1(α::Integer, β::Integer, γ::Integer, δ::Integer; check_args::Bool=true) = Stable1(float(α), float(β), float(γ), float(δ); check_args=check_args)
Stable1(α::Real, β::Real; check_args::Bool=true) = Stable1(α, β, one(α), zero(α); check_args=check_args)
Stable1(α::Real; check_args::Bool=true) = Stable1(α, zero(α); check_args=check_args)


@distr_support Stable1 (d.α < 1 && d.β == 1 ? d.μ1 : -Inf) (d.α < 1 && d.β == -1 ? d.μ1 : Inf)

#### Conversions

convert(::Type{Stable1{T}}, α::S, β::S, γ::S, μ::S) where {T <: Real, S <: Real} = Stable1(T(α), T(β), T(γ), T(δ))
Base.convert(::Type{Stable1{T}}, d::Stable1) where {T<:Real} = Stable1{T}(T(d.α), T(d.β), T(d.γ), T(d.δ))
Base.convert(::Type{Stable1{T}}, d::Stable1{T}) where {T<:Real} = d

#### Parameters

shape(d::Stable1) = (d.α, d.β)
location(d::Stable1) = d.γ
scale(d::Stable1) = d.δ
params(d::Stable1) = (d.α, d.β, d.γ, d.δ)
partype(::Stable1{T}) where {T} = T

parametrization(d::Stable1)= 1
params1(d::Stable1) = (d.α, d.β, d.σ1, d.μ1)

#### Statistics

mean(d::Stable1{T}) where T = d.α > one(T) ? d.μ1 : T(NaN)
var(d::Stable1{T}) where T = d.α == 2one(T) ? 2d.σ1^2 : T(Inf)
skewness(d::Stable1{T}) where T = d.α == 2one(T) ? T(0.0) : T(NaN)
kurtosis(d::Stable1{T}) where T = d.α == 2one(T) ? T(0.0) : T(NaN)

#### Evaluation

# ok for parametrization 1, to be checked for parametrization 0 and 2
function cf(d::Stable1{T}, t::Real) where T
    ds=Stable(d.α1, d.β1, d.σ1, d.μ1)
    return cf(ds,x)
end

# ok for parametrization 1, to be checked for parametrization 0 and 2
function mgf(d::Stable1{T}, t::Real) where T
    ds=Stable(d.α1, d.β1, d.σ1, d.μ1)
    return mgf(ds,x)
end

function pdf(d::Stable1{T}, x::Real) where T
    ds=Stable(d.α1, d.β1, d.σ1, d.μ1)
    return pdf(ds,x)
end

function logpdf(d::Stable1{T}, x::Real) where T
    ds=Stable(d.α1, d.β1, d.σ1, d.μ1)
    return logpdf(ds,x)
end

function cdf(d::Stable1{T}, x::Real) where T
    ds=Stable(d.α1, d.β1, d.σ1, d.μ1)
    return cdf(ds,x)
end

# """
#    mode : compute mode for for Stable distribution(α1,β1,γ1,δ1) with parametrization 1
#    see Nolan's chapter 3.2, p.77
# """
function mode(d::Stable1{T}) where T
    α, β, γ, δ = (d.α, d.β, d.γ, d.δ)
    d=0.0
    if !(α ≈ 1.0)
        d=δ+β*γ*tanpi(α/2)
    else
        d=δ+β*2/π*γ*log(γ)
    end
    return γ*_mm(α,β)+d
end

quantile(d::Stable1, p::Real)   = Stable.quantile_newton(StablefromStable1(d), p, mode(d), 1e-6)
cquantile(d::Stable1, p::Real)  = Stable.cquantile_newton(StablefromStable1(d), p, mode(d), 1e-6)
invlogcdf(d::Stable1, p::Real)  = Stable.invlogcdf_newton(StablefromStable1(d), p, mode(d), 1e-6)
invlogccdf(d::Stable1, p::Real) = Stable.invlogccdf_newton(StablefromStable1(d), p, mode(d), 1e-6)

#### Affine transformations

Base.:+(d::Stable1, a::Real) = Stable1(d.α, d.β, d.γ, d.δ + a)
Base.:*(c::Real, d::Stable1{T}) where T = Stable1fromStable(c*StablefromStable1(d))

#### Sampling

# A. Weron, R. Weron "Computer simulation of Lévy α-stable variables and processes", Springer 1995, doi: 10.1007/3-540-60188-0_6
function rand(rng::AbstractRNG, d::Stable1{T}) where T
    (α, β, σ, μ) = (d.α1, d.β1, d.σ1, d.μ1)
    return rand(rng, Stable(d.α1, d.β1, d.σ1, d.μ1))
end


#### Fit model

"""
    fit_quantile(::Type{<:Stable1}, x::AbstractArray{<:Real})

Compute the estimate of the [`Stable`](@ref) distribution with McCulloch's quantile method.
    Result is tuple with distribution's type-0 parameters.
"""
function fit_quantile(::Type{<:Stable1}, x::AbstractArray{<:Real})
    α₀, β₀, σ₀, δ₀ = fit_quantile(Stable, x)
    return (α₀, β₀, σ₀, δ₀)
end

"""
    fit(::Type{<:Stable1}, x::AbstractArray{<:Real})
    ECF method, see Nolan ch. 4.
    quantile_fit is used for a robust initial estimate
"""
function fit(::Type{<:Stable1}, x::AbstractArray{<:Real})
    dStablefitted=fit(Stable, x)
    (α, β, γ, δ) = _paramconvstabledists(1,2,params(dStablefitted)...)
    return Stable1(α, β, γ, δ)
end
