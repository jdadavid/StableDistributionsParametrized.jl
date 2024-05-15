# extension to StableDistributions to using Nolan's parametrizations (0,1,2)
# knowing that StableDistributions implicit parametrization is 1, whereas Nolan uses 0 by default in his book

# methods mainly from John P. Nolan, "Univariate Stable Distributions", Springer 2020

"""
    Stable2(α, β, γ, δ)

The *Stable distribution* with stability index 0 < `α` ≤ 2, skewness parameter -1 ≤ `β` ≤ 1, scale 0 < `γ` and location `δ`
with Nolan'sparametrization 2

From from John P. Nolan, "Univariate Stable Distributions", Springer 2020, chapter 1.3 (p.7)

    p=0  is best for numerical works and statisticsinference; the characteristic function cf is continuous in all parameters

    p=1  is best for working with characteristic function and nice algebraic properties
         but location of the mode is unbounded near α=1 (tends to +Inf if α<1 and to -Inf if α>1)

    p=2  a.k.a mode-centered, is best for users
         location and scale parameters agree with standard Normal's when α=2, and with Cauchy's when α=1

```julia

Stable2(α)            # standard symmetric α-Stable distribution equivalent to Stablep(α, 0, 1, 0)
Stable2(α, β)         # standard α-Stable distribution with skewness parameter β equivalent to S(α, β, 1, 0)
Stable2(α, β, γ, δ)   # α-Stable distribution with skewness parameter β, scale γ and location δ

params(d)             # Get the parameters, i.e. (α, β, γ, δ)
shape(d)              # Get the shape, i.e. (α, β)
location(d)           # Get the location, i.e. δ
scale(d)              # Get the scale, i.e. γ
minimum(d)            # Get the lower bound
maximum(d)            # Get the upper bound
```

External links

* [Stable distribution on Wikipedia](https://en.wikipedia.org/wiki/Stable_distribution)

"""
struct Stable2{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    β::T
    γ::T
    δ::T
    α1::T
    β1::T
    σ1::T
    μ1::T
    function Stable2{T}(α, β, γ, δ; p, check_args::Bool=true) where {T}
        @check_args Stable (α, zero(T) < α <= 2one(T)) (β, -one(T) <= β <= one(T)) (γ, zero(T) < γ) ((α, β), α != 2 || β == 0 )
        (α1,β1,σ1,μ1) = _paramconvstabledists(2,1,α,β,γ,δ)
        new{T}(α, β, γ, δ, α1, β1, σ1, μ1)
    end
end

"""
    Stable1(α, β, γ, δ)

The *Stable distribution* with stability index 0 < `α` ≤ 2, skewness parameter -1 ≤ `β` ≤ 1, scale 0 < `γ` and location `δ`
with Nolan'sparametrization 1
Same as Stable

From from John P. Nolan, "Univariate Stable Distributions", Springer 2020, chapter 1.3 (p.7)

    p=0  is best for numerical works and statisticsinference; the characteristic function cf is continuous in all parameters

    p=1  is best for working with characteristic function and nice algebraic properties
         but location of the mode is unbounded near α=1 (tends to +Inf if α<1 and to -Inf if α>1)

    p=2  a.k.a mode-centered, is best for users
         location and scale parameters agree with standard Normal's when α=2, and with Cauchy's when α=1

```julia

Stable1(α)            # standard symmetric α-Stable distribution equivalent to Stablep(α, 0, 1, 0)
Stable1(α, β)         # standard α-Stable distribution with skewness parameter β equivalent to S(α, β, 1, 0)
Stable1(α, β, γ, δ)   # α-Stable distribution with skewness parameter β, scale γ and location δ

params(d)             # Get the parameters, i.e. (α, β, γ, δ)
shape(d)              # Get the shape, i.e. (α, β)
location(d)           # Get the location, i.e. δ
scale(d)              # Get the scale, i.e. γ
minimum(d)            # Get the lower bound
maximum(d)            # Get the upper bound
```

External links

* [Stable distribution on Wikipedia](https://en.wikipedia.org/wiki/Stable_distribution)

"""
struct Stable2{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    β::T
    γ::T
    δ::T
    α1::T
    β1::T
    σ1::T
    μ1::T
    function Stable2{T}(α, β, γ, δ; p, check_args::Bool=true) where {T}
        @check_args Stable (α, zero(T) < α <= 2one(T)) (β, -one(T) <= β <= one(T)) (γ, zero(T) < γ) ((α, β), α != 2 || β == 0 )
        # (α1,β1,σ1,μ1)=_paramconvstabledists(1,1,α,β,γ,δ)
        (α1,β1,σ1,μ1) = (α1,β1,σ1,μ1)
        new{T}(α, β, γ, δ, α1, β1, σ1, μ1)
    end
end

"""
    Stable0(α, β, γ, δ)

The *Stable distribution* with stability index 0 < `α` ≤ 2, skewness parameter -1 ≤ `β` ≤ 1, scale 0 < `γ` and location `δ`
with Nolan'sparametrization 0

From from John P. Nolan, "Univariate Stable Distributions", Springer 2020, chapter 1.3 (p.7)

    p=0  is best for numerical works and statisticsinference; the characteristic function cf is continuous in all parameters

    p=1  is best for working with characteristic function and nice algebraic properties
         but location of the mode is unbounded near α=1 (tends to +Inf if α<1 and to -Inf if α>1)

    p=2  a.k.a mode-centered, is best for users
         location and scale parameters agree with standard Normal's when α=2, and with Cauchy's when α=1

```julia

Stable0(α)            # standard symmetric α-Stable distribution equivalent to Stablep(α, 0, 1, 0)
Stable0(α, β)         # standard α-Stable distribution with skewness parameter β equivalent to S(α, β, 1, 0)
Stable0(α, β, γ, δ)   # α-Stable distribution with skewness parameter β, scale γ and location δ

params(d)             # Get the parameters, i.e. (α, β, γ, δ)
shape(d)              # Get the shape, i.e. (α, β)
location(d)           # Get the location, i.e. δ
scale(d)              # Get the scale, i.e. γ
minimum(d)            # Get the lower bound
maximum(d)            # Get the upper bound
```

External links

* [Stable distribution on Wikipedia](https://en.wikipedia.org/wiki/Stable_distribution)

"""
struct Stable0{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    β::T
    γ::T
    δ::T
    α1::T
    β1::T
    σ1::T
    μ1::T
    function Stable0{T}(α, β, γ, δ; p, check_args::Bool=true) where {T}
        @check_args Stable (α, zero(T) < α <= 2one(T)) (β, -one(T) <= β <= one(T)) (γ, zero(T) < γ) ((α, β), α != 2 || β == 0 )
        (α1,β1,σ1,μ1)=_paramconvstabledists(0,1,α,β,γ,δ)
        new{T}(α, β, γ, δ, α1, β1, σ1, μ1)
    end
end

# """
#     _paramconvstabledists(pfrom,pto,α,β,γ,δ) : convert parametrization "pfrom" to parametrization "pto"
#       see table 3.1 of Chapter 3.5 (p.97) of Univariate Stable Distributions - John P. Nolan (2020)
# """
# function _paramconvstabledists(pfrom,pto,α,β,γ,δ)
#     if pfrom == pto ; return α,β,γ,δ; end
#     if !(pfrom in [0,1,2])
#         error("stabledistributionsx : stabledists_paramconv : pfrom=$pfrom not in [0,1,2]")
#     end
#     if !(pto in [0,1,2])
#         error("stabledistributionsx : stabledists_paramconv : pto=$pto not in [0,1,2]")
#     end
#     (a,b,c,d)=(α,β,γ,δ)
#
#     if     pfrom == 0
#         if     pto == 1
#             if !(α ≈ 1.0)
#                 d=δ-β*γ*tanpi(α/2)
#             else
#                 d=δ-β*2/π*γ*log(γ)
#             end
#         elseif pto == 2
#             c=γ/(α^(-1.0/α))
#             d=δ+γ*_mm(α,β)
#         end
#     elseif pfrom == 1
#         if     pto == 0
#             if !(α ≈ 1.0)
#                 d=δ+β*γ*tanpi(α/2)
#             else
#                 d=δ+β*2/π*γ*log(γ)
#             end
#         elseif pto == 2
#             c=γ/(α^(-1.0/α))
#             if !(α ≈ 1.0)
#                 d=δ+γ*(β*tanpi(α/2)+_mm(α,β))
#             else
#                 d=δ+c*(β*2/π*log(c)+_mm(α,β))
#             end
#         end
#     elseif pfrom == 2
#         c=α^(-1.0/α)*γ
#         if     pto == 0
#             d=δ-c*mm(α,β)
#         elseif pto == 1
#             if !(α ≈ 1.0)
#                 d=δ-c*(β*tanpi(α/2)+_mm(α,β))
#             else
#                 d=δ-γ*(β*2/π*log(γ)+_mm(α,β))
#             end
#         end
#     end
#
#     return (a,b,c,d)
# end

StablefromStable2(dp::Stablep)   = Stable(dp.α1, dp.β1, dp.σ1, dp.μ1)
Stable2fromStable(d::Stable)     = Stable2(_paramconvstabledists(1,2,params(d)...)...)

StablefromStable1(dp::Stablep)   = Stable(dp.α1, dp.β1, dp.σ1, dp.μ1)
Stable1fromStable(d::Stable)     = Stable1(_paramconvstabledists(1,2,params(d)...)...)

StablefromStable0(dp::Stablep)   = Stable(dp.α1, dp.β1, dp.σ1, dp.μ1)
Stable0fromStable(d::Stable)     = Stable2(_paramconvstabledists(1,0,params(d)...)...)


Stable2(α, β)       = Stablep2(α, β  , 1.0,0.0)
Stable2(α,)         = Stablep2(α, 0.0, 1.0,0.0)

Stable1(α, β)       = Stablep1(α, β  , 1.0,0.0)
Stable1(α,)         = Stablep1(α, 0.0, 1.0,0.0)


Stable0(α, β)       = Stablep0(α, β  , 1.0,0.0)
Stable0(α,)         = Stablep0(α, 0.0, 1.0,0.0)

Stablep2fromStable(d::Stable) = Stablep(_paramconvstabledists(1,p,params(d)...)...;p=2)
Stablep1fromStable(d::Stable) = Stablep(_paramconvstabledists(1,p,params(d)...)...;p=1)
Stablep0fromStable(d::Stable) = Stablep(_paramconvstabledists(1,p,params(d)...)...;p=0)

StablefromStablep2(dp::Stablep)   = Stable(dp.α1, dp.β1, dp.σ1, dp.μ1)
StablefromStablep1(dp::Stablep)   = Stable(dp.α1, dp.β1, dp.σ1, dp.μ1)
StablefromStablep0(dp::Stablep)   = Stable(dp.α1, dp.β1, dp.σ1, dp.μ1)


Stablep(α::T, β::T, γ::T, δ::T; p=1, check_args::Bool=true) where {T <: Real} = Stablep{T}(α, β, γ, δ; p=p, check_args=check_args)

Stablep(α::Real, β::Real, γ::Real, δ::Real; p=1, check_args::Bool=true) = Stablep(promote(α, β, γ, δ)...; p=p, check_args=check_args)
Stablep(α::Integer, β::Integer, γ::Integer, δ::Integer; check_args::Bool=true) = Stablep(float(α), float(β), float(γ), float(δ); p=p, check_args=check_args)
Stablep(α::Real, β::Real; p=1, check_args::Bool=true) = Stablep(α, β, one(α), zero(α); p=p, check_args=check_args)
Stablep(α::Real; p=1, check_args::Bool=true) = Stable(α, zero(α); p=p, check_args=check_args)


@distr_support Stablep (d.α < 1 && d.β == 1 ? d.μ1 : -Inf) (d.α < 1 && d.β == -1 ? d.μ1 : Inf)

#### Conversions

convert(::Type{Stable2{T}}, α::S, β::S, γ::S, μ::S) where {T <: Real, S <: Real} = Stable2(T(α), T(β), T(γ), T(δ))
Base.convert(::Type{Stable2{T}}, d::Stable2) where {T<:Real} = Stable2{T}(T(d.α), T(d.β), T(d.γ), T(d.δ))
Base.convert(::Type{Stable2{T}}, d::Stable2{T}) where {T<:Real} = d

#### Parameters

shape(d::Stable2) = (d.α, d.β)
location(d::Stable2) = d.μ1
scale(d::Stable2) = d.σ1
params(d::Stable2) = (d.α, d.β, d.γ, d.δ)
partype(::Stable2{T}) where {T} = T

parametrization(d::Stable2)= 2
params1(d::Stable2) = (d.α, d.β, d.σ1, d.μ1)

#### Statistics

mean(d::Stable2{T}) where T = d.α > one(T) ? d.μ1 : T(NaN)
var(d::Stable2{T}) where T = d.α == 2one(T) ? 2d.σ1^2 : T(Inf)
skewness(d::Stable2{T}) where T = d.α == 2one(T) ? T(0.0) : T(NaN)
kurtosis(d::Stable2{T}) where T = d.α == 2one(T) ? T(0.0) : T(NaN)

#### Evaluation

# ok for parametrization 1, to be checked for parametrization 0 and 2
function cf(d::Stable2{T}, t::Real) where T
    ds=Stable(d.α1, d.β1, d.σ1, d.μ1)
    return cf(ds,x)
#     # α, β, σ, μ =  params(d)
#     α, β, σ, μ = (d.α1, d.β1, d.σ1, d.μ1)
#     if α == one(T)
#         exp(im*t*μ - abs(σ*t) * (1 + im*β*2/π*sign(t)*log(abs(t))))
#     else
#         exp(im*t*μ - abs(σ*t)^α * (1 - im*β*sign(t)*tan(α*π/2)))
#     end
end

# ok for parametrization 1, to be checked for parametrization 0 and 2
function mgf(d::Stable2{T}, t::Real) where T
    ds=Stable(d.α1, d.β1, d.σ1, d.μ1)
    return mgf(ds,x)
#     if d.α1 == 2one(T)
#         mgf(Normal(d.μ, √2d.σ), t)
#     else
#         T(Inf)
#     end
end

function pdf(d::Stable2{T}, x::Real) where T
    ds=Stable(d.α1, d.β1, d.σ1, d.μ1)
    return pdf(ds,x)
end

function logpdf(d::Stable2{T}, x::Real) where T
    ds=Stable(d.α1, d.β1, d.σ1, d.μ1)
    return logpdf(ds,x)
end

function cdf(d::Stable2{T}, x::Real) where T
    ds=Stable(d.α1, d.β1, d.σ1, d.μ1)
    return cdf(ds,x)
end

# """
#     _mm(α,β) : compute mode for Stable distribution(α,β) with parametrization 0
#       see Appendix C of Univariate Stable Distributions - John P. Nolan (2020) for data
#       function is using  linear interpolation (shameless copied from StableDistributions src) between Nolan's tabulated data
# """
# function _mm(α,β)
#     if !(0.01 <= α <= 2.0)
#         error("Stablep : _mm : α=$α not inside ]0.01,2.0]")
#     end
#     if !(-1.0 <= β <= 1.0)
#         error("Stablep : _mm : β=$β not inside [-1.0,1.0]")
#     end
#     alpha=α
#     beta=β
#     if β<0 ; betaneg=true; beta=-β; else; betaneg=false; end
#
#     xs=0.0:0.05:2.0 # alpha range
#     ys=0.0:0.1:1.0  # beta  range
#
# # from Appendix C of Univariate Stable Distributions - John P. Nolan (2020)
# # Line for α=0 contains invalid values, but is used for computing index for interpolation
#  #    β =   0.0      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0   #  α
#     A  = [0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000 # 0.00 INVALID?
#           0.00000 -0.00787 -0.01574 -0.02361 -0.03148 -0.03935 -0.04722 -0.05509 -0.06296 -0.07083 -0.07870 # 0.05
#           0.00000 -0.01584 -0.03168 -0.04752 -0.06335 -0.07919 -0.09503 -0.11087 -0.12671 -0.14255 -0.15338 # 0.10
#           0.00000 -0.02401 -0.04802 -0.07202 -0.09603 -0.12004 -0.14405 -0.16806 -0.19206 -0.21607 -0.23508 # 0.15
#           0.00000 -0.03249 -0.06498 -0.09748 -0.12997 -0.16246 -0.19495 -0.22744 -0.25994 -0.28743 -0.31992 # 0.20
#           0.00000 -0.04142 -0.08284 -0.12426 -0.16569 -0.20711 -0.24353 -0.28495 -0.32637 -0.36779 -0.40921 # 0.25
#           0.00000 -0.05095 -0.10191 -0.14786 -0.19881 -0.24976 -0.29986 -0.34879 -0.39731 -0.44539 -0.49299 # 0.30
#           0.00000 -0.06128 -0.11756 -0.17705 -0.23479 -0.29177 -0.34792 -0.40319 -0.45751 -0.51082 -0.56307 # 0.35
#           0.00000 -0.06765 -0.13452 -0.20004 -0.26418 -0.32686 -0.38800 -0.44754 -0.50540 -0.56153 -0.61587 # 0.40
#           0.00000 -0.07454 -0.14754 -0.21852 -0.28730 -0.35380 -0.41794 -0.47967 -0.53894 -0.59572 -0.64995 # 0.45
#           0.00000 -0.07955 -0.15704 -0.23174 -0.30340 -0.37193 -0.43729 -0.49945 -0.55841 -0.61415 -0.66667 # 0.50
#           0.00000 -0.08267 -0.16288 -0.23963 -0.31258 -0.38163 -0.44678 -0.50804 -0.56544 -0.61902 -0.66880 # 0.55
#           0.00000 -0.08399 -0.16527 -0.24259 -0.31549 -0.38388 -0.44776 -0.50719 -0.56223 -0.61297 -0.65948 # 0.60
#           0.00000 -0.08372 -0.16463 -0.24127 -0.31307 -0.37990 -0.44177 -0.49876 -0.55099 -0.59857 -0.64162 # 0.65
#           0.00000 -0.08213 -0.16147 -0.23642 -0.30630 -0.37092 -0.43030 -0.48452 -0.53372 -0.57806 -0.61768 # 0.70
#           0.00000 -0.07949 -0.15631 -0.22878 -0.29613 -0.35810 -0.41468 -0.46597 -0.51211 -0.55327 -0.58963 # 0.75
#           0.00000 -0.07606 -0.14964 -0.21904 -0.28341 -0.34243 -0.39605 -0.44436 -0.48750 -0.52565 -0.55899 # 0.80
#           0.00000 -0.07208 -0.14188 -0.20777 -0.26886 -0.32474 -0.37534 -0.42070 -0.46097 -0.49631 -0.52692 # 0.85
#           0.00000 -0.06772 -0.13339 -0.19549 -0.25308 -0.30572 -0.35328 -0.39577 -0.43332 -0.46609 -0.49424 # 0.90
#           0.00000 -0.06314 -0.12448 -0.18259 -0.23656 -0.28592 -0.33047 -0.37020 -0.40519 -0.43559 -0.46156 # 0.95
#           0.00000 -0.05847 -0.11537 -0.16940 -0.21970 -0.26577 -0.30736 -0.34443 -0.37704 -0.40528 -0.42931 # 1.00
#           0.00000 -0.05380 -0.10625 -0.15619 -0.20281 -0.24561 -0.28431 -0.31884 -0.34920 -0.37547 -0.39778 # 1.05
#           0.00000 -0.04921 -0.09725 -0.14315 -0.18613 -0.22571 -0.26160 -0.29367 -0.32192 -0.34639 -0.36717 # 1.10
#           0.00000 -0.04474 -0.08850 -0.13042 -0.16983 -0.20626 -0.23941 -0.26914 -0.29539 -0.31818 -0.33758 # 1.15
#           0.00000 -0.04043 -0.08005 -0.11812 -0.15406 -0.18742 -0.21791 -0.24537 -0.26972 -0.29094 -0.30909 # 1.20
#           0.00000 -0.03632 -0.07197 -0.10633 -0.13891 -0.16929 -0.19720 -0.22247 -0.24499 -0.26474 -0.28173 # 1.25
#           0.00000 -0.03242 -0.06429 -0.09511 -0.12444 -0.15195 -0.17736 -0.20051 -0.22127 -0.23961 -0.25550 # 1.30
#           0.00000 -0.02873 -0.05703 -0.08448 -0.11072 -0.13545 -0.15845 -0.17953 -0.19859 -0.21555 -0.23040 # 1.35
#           0.00000 -0.02527 -0.05020 -0.07446 -0.09775 -0.11983 -0.14049 -0.15957 -0.17696 -0.19258 -0.20639 # 1.40
#           0.00000 -0.02204 -0.04381 -0.06507 -0.08556 -0.10510 -0.12350 -0.14063 -0.15639 -0.17068 -0.18347 # 1.45
#           0.00000 -0.01903 -0.03786 -0.05629 -0.07415 -0.09126 -0.10750 -0.12273 -0.13687 -0.14984 -0.16159 # 1.50
#           0.00000 -0.01624 -0.03233 -0.04813 -0.06350 -0.07832 -0.09248 -0.10587 -0.11842 -0.13006 -0.14073 # 1.55
#           0.00000 -0.01366 -0.02722 -0.04058 -0.05362 -0.06627 -0.07843 -0.09003 -0.10101 -0.11131 -0.12088 # 1.60
#           0.00000 -0.01130 -0.02252 -0.03361 -0.04448 -0.05508 -0.06534 -0.07522 -0.08465 -0.09360 -0.10202 # 1.65
#           0.00000 -0.00913 -0.01822 -0.02721 -0.03607 -0.04475 -0.05321 -0.06142 -0.06933 -0.07692 -0.08415 # 1.70
#           0.00000 -0.00716 -0.01429 -0.02137 -0.02837 -0.03526 -0.04202 -0.04863 -0.05505 -0.06128 -0.06729 # 1.75
#           0.00000 -0.00538 -0.01074 -0.01607 -0.02136 -0.02660 -0.03176 -0.03684 -0.04183 -0.04670 -0.05146 # 1.80
#           0.00000 -0.00377 -0.00754 -0.01130 -0.01503 -0.01874 -0.02243 -0.02607 -0.02967 -0.03323 -0.03673 # 1.85
#           0.00000 -0.00235 -0.00469 -0.00703 -0.00937 -0.01170 -0.01401 -0.01632 -0.01861 -0.02089 -0.02315 # 1.90
#           0.00000 -0.00109 -0.00218 -0.00327 -0.00436 -0.00545 -0.00653 -0.00762 -0.00870 -0.00978 -0.01086 # 1.95
#           0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000]# 2.00
#
#     valmm=0.0
#
#     function interpolate2D(x, y, xs, ys, M)
#         lx, ly = length(xs), length(ys)
#         i, j  = findlast(xs .<= x), findlast(ys .<= y)
#         i == lx && j < ly && return M[end, j]
#         j == ly && i < lx && return M[i, end]
#         i == lx && j == ly && return M[end,end]
#         x₁, x₂, y₁, y₂ = xs[i], xs[i+1], ys[j], ys[j+1]
#         m₁₁, m₁₂, m₂₁, m₂₂ = M[i, j], M[i, j+1], M[i+1, j], M[i+1, j+1]
#         xmy = (x₂-x) * (m₁₁*(y₂-y) + m₁₂*(y-y₁)) + (x-x₁) * (m₂₁*(y₂-y) + m₂₂*(y-y₁))
#         return xmy / ( (x₂-x₁)*(y₂-y₁) )
#     end
#
#     valmm=interpolate2D(alpha,beta,xs,ys,A)
#
#     if betaneg; valmm=-valmm; end
#     return valmm
# end

"""
   mode : compute mode for for Stable distribution(α1,β1,γ1,δ1) with parametrization 1
   see Nolan's chapter 3.2, p.77
"""
function mode(d::Stablep{T}) where T
    α, β, γ, δ = (d.α, d.β, d.γ, d.δ)
    d=0.0
    if !(α ≈ 1.0)
        d=δ+β*γ*tanpi(α/2)
    else
        d=δ+β*2/π*γ*log(γ)
    end
    return γ*_mm(α,β)+d
end

quantile(d::Stablep, p::Real)   = Stable.quantile_newton(StablefromStablep(d), p, mode(d), 1e-6)
cquantile(d::Stablep, p::Real)  = Stable.cquantile_newton(StablefromStablep(d), p, mode(d), 1e-6)
invlogcdf(d::Stablep, p::Real)  = Stable.invlogcdf_newton(StablefromStablep(d), p, mode(d), 1e-6)
invlogccdf(d::Stablep, p::Real) = Stable.invlogccdf_newton(StablefromStablep(d), p, mode(d), 1e-6)

#### Affine transformations

Base.:+(d::Stablep, a::Real) = Stablep(d.α, d.β, d.γ, d.δ + a)
# Base.:*(c::Real, d::Stablep{T}) where T =
#     if d.α == one(T)
#         Stablep(d.α, sign(c)*d.β, abs(c)*d.γ, c*d.δ - 2/π*d.β*d.γ*c*log(abs(c)))
#     else
#         Stablep(d.α, sign(c)*d.β, abs(c)*d.γ, c*d.δ)
#     end
Base.:*(c::Real, d::Stablep{T}) where T = StablepfromStable(c*StablefromStablep(d);p=d.p)

#### Sampling

# A. Weron, R. Weron "Computer simulation of Lévy α-stable variables and processes", Springer 1995, doi: 10.1007/3-540-60188-0_6
function rand(rng::AbstractRNG, d::Stablep{T}) where T
    (α, β, σ, μ) = (d.α1, d.β1, d.σ1, d.μ1)
    return rand(rng, Stable(d.α1, d.β1, d.σ1, d.μ1))
end


#### Fit model

"""
    fit_quantile(::Type{<:Stable}, x::AbstractArray{<:Real})

Compute the estimate of the [`Stable`](@ref) distribution with McCulloch's quantile method.
    Result is tuple with distribution's type-0 parameters.
"""
function fit_quantile(::Type{<:Stablep}, x::AbstractArray{<:Real})
    α₀, β₀, σ₀, δ₀ = fit_quantile(Stable, x)
    return (α₀, β₀, σ₀, δ₀)
end

"""
    fit(::Type{<:Stablep}, x::AbstractArray{<:Real};p=1)
    ECF method, see Nolan ch. 4.
    quantie_fit is used for a robust initial estimate
"""
function fit(::Type{<:Stablep}, x::AbstractArray{<:Real};p=1)
    dStablefitted=fit(Stable, x)
    (α, β, γ, δ) = _paramconvstabledists(1,p,dStablefitted)
    return Stablep(α, β, γ, δ; p=p)
end
