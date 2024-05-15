using Distributions
using StableDistributions
if !isdefined(Main,:StableDistributionsp); include("StableDistributionsp.jl"); using ..StableDistributionsp; end
using StatsBase
using StatsPlots
# using Interpolations
# using Optim
# using RCall
# R"library(stabledist)"

"""
    plotnolan_fig1_2 : plot Nolan' fig 1.2, "Stable Densities S(α,0.5;0), parametrization 0, α=0.5,0.75,1,1.25,1.5"
      see figure 1.2 of Chapter 1 (p.8) of Univariate Stable Distributions - John P. Nolan (2020)
"""
function plot_nolan_fig1_2()
    p=plot(title="Stable Densities S(α,0.5) parametrization 0",xlims=(-5.0,5.0),xlabel="α")
    β=0.5
    for α in (0.5,0.75,1.0,1.25,1.5)
        # (a,b,c,d)=_paramconvstabledists(0,1,α,β,1.0,0.0)
        # d=Stable(a,b,c,d)
        d=Stablep(α,β,1.0,0.0;p=0)
        p=plot!(p,x->pdf(d,x),label="α=$α")
    end
    display(p)
    return nothing
end

"""
    plotnolan_fig1_3 : plot Nolan' fig 1.3, "Stable Densities S(α,0.5;1), parametrization 1, α=0.5,0.75,1,1.25,1.5"
      see figure 1.3 of Chapter 1 (p.8) of Univariate Stable Distributions - John P. Nolan (2020)
"""
function plot_nolan_fig1_3()
    p=plot(title="Stable Densities S(α,0.5) parametrization 1",xlims=(-5.0,5.0),xlabel="α")
    β=0.5
    for α in (0.5,0.75,1.0,1.25,1.5)
        # d=Stable(α,β)
        d=Stablep(α,β,1.0,0.0;p=1)
        p=plot!(p,x->pdf(d,x),label="α=$α")
    end
    display(p)
    return nothing
end

"""
    plotnolan_fig1_4 : plot Nolan' fig 1.4, "Stable Densities S(α,0.5;2), parametrization 2, α=0.25,0.75,1,1.25,1.5"
      see figure 1.4 of Chapter 1 (p.8) of Univariate Stable Distributions - John P. Nolan (2020)
"""
function plot_nolan_fig1_4()
    p=plot(title="Stable Densities S(α,0.5,1) parametrization 2",xlims=(-5.0,5.0),xlabel="α")
    β=0.5
    for α in (0.5,0.75,1.0,1.25,1.5)
        # (a,b,c,d)=_paramconvstabledists(2,1,α,β,1.0,0.0)
        # d=Stable(a,b,c,d)
        d=Stablep(α,β,1.0,0.0;p=2)
        p=plot!(p,x->pdf(d,x),label="α=$α")
    end
    display(p)
    return nothing
end

function pb0old()
    (a,b,c,d)=_paramconvstabledists(0,1,0.75,0.5,1.0,0.0)
    # (0.75, 0.5, 1.0, -1.2071067811865475)
    d=Stable(a,b,c,d)
    plot(x->pdf(d,x))
end
function pb2old()
    (a,b,c,d)=_paramconvstabledists(2,1,0.75,0.5,1.0,0.0)
    # (0.75, 0.5, 1.467523221730945, -1.2459371667983017)
    d=Stable(a,b,c,d)
    plot(x->pdf(d,x))
end

function pb0()
    #(a,b,c,d)=_paramconvstabledists(0,1,0.75,0.5,1.0,0.0)
    # (0.75, 0.5, 1.0, -1.2071067811865475)
    d=Stablep(0.75,0.5,1.0,0.0;p=0)
    plot(x->pdf(d,x))
end
function pb2()
    #(a,b,c,d)=_paramconvstabledists(2,1,0.75,0.5,1.0,0.0)
    # (0.75, 0.5, 1.467523221730945, -1.2459371667983017)
    d=Stablep(0.75,0.5,1.0,0.0;p=2)
    plot(x->pdf(d,x))
end

nothing
