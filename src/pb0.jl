using Distributions
using StableDistributions
if !isdefined(Main,:StableDistributionsp); include("StableDistributionsp.jl"); using ..StableDistributionsp; end
using StatsBase
using StatsPlots
# using Interpolations
# using Optim
# using RCall
# R"library(stabledist)"

function pb0()
    #(a,b,c,d)=_paramconvstabledists(0,1,0.75,0.5,1.0,0.0)
    # (0.75, 0.5, 1.0, -1.2071067811865475)
    d=Stablep(0.75,0.5,1.0,0.0;p=0)
    plot(x->pdf(d,x))
end
function pb0re()
    #(a,b,c,d)=_paramconvstabledists(0,1,0.75,0.5,1.0,0.0)
    # (0.75, 0.5, 1.0, -1.2071067811865475)
    d=Stablep(0.75,0.5,1.0,0.0;p=0)
    dpdf(x)=pdf(d,x)
    p=plot(title="pb0re")
    # p=plot!(p,dpdf,label="dpdf")
    xr=-5.0:0.01:5.0
    xr=-2.0:0.01:-1.0
    xr=-1.25:0.01:-1.1
    xr=-1.23:0.001:-1.17
    xr=-1.210:0.0001:-1.194
    yd=dpdf.(xr)
    p=plot!(p,xr,yd,label="xryr")
    display(p)
end
function pb0re1()
    d=Stable(0.75,0.5,1.0,0.0)
    dpdf(x)=pdf(d,x)
    p=plot(title="pb0re1")
    xr=-0.010:0.00001:0.02
    yd=dpdf.(xr)
    p=plot!(p,xr,yd,label="xryr")
    display(p)
end
function pb0mwe()
    d=Stable(0.75,0.5,1.0,0.0)
    x=0.005
    y=pdf(d,x)
    println("x=",x," y=",y)
end
nothing

# pb0re()
pb0mwe()
pb0re1()
