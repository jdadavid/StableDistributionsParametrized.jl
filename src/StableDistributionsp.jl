module StableDistributionsp

using Random
import Random: rand, AbstractRNG
using Distributions
using StableDistributions
# using SpecialFunctions
# using QuadGK: quadgk
# using StatsFuns: logexpm1, log1mexp
import StableDistributions: @check_args, @distr_support,
    params, shape, location, scale, support, minimum, maximum,
    convert, convolve, +, *,
    partype, mean, var, skewness, kurtosis,
    mgf, cf, pdf, logpdf, cdf, fit,
    quantile, cquantile, invlogcdf, invlogccdf


export Stablep,                                # starting specific exports for Stablep
    Stablep2,Stablep1,Stablep0,
    Stable2,Stable1,Stable0,                   # shorthands for Stablep(...;p=0|1|2)
    _paramconvstabledists,
    StablefromStablep, StablepfromStable,
    StablefromStablep2, Stablep2fromStable,
    StablefromStablep1, Stablep1fromStable,
    StablefromStablep0, Stablep0fromStable,
    StablefromStable2, Stable2fromStable,
    StablefromStable1, Stable1fromStable,
    StablefromStable0, Stable0fromStable,
    StablepfromStable2, StablepfromStable1, StablepfromStable0,
    Stable2fromStablep, Stable1fromStablep, Stable0fromStablep,
    mode,                                      # note : mode for Stable is defined in include("modeforStable.jl")
    paramsp,parametrization, paramsp1,
    Stable,                                    # starting similar exports as Stable
    rand,
    params, shape, location, scale, support, minimum, maximum,
    convert, convolve, +, *,
    partype, mean, var, skewness, kurtosis,
    mgf, cf, pdf, logpdf, cdf, mgf, fit, fit_quantile,
    quantile, cquantile, invlogcdf, invlogccdf



include("stablep.jl")
include("modeforStable.jl") # separate file to note extension to Stable
include("stable2.jl")
include("stable1.jl")
include("stable0.jl")

end
