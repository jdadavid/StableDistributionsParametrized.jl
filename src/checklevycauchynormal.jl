using Distributions
using StableDistributions
if !isdefined(Main,:StableDistributionsp); include("StableDistributionsp.jl"); using ..StableDistributionsp; end
using StatsBase
using StatsPlots
# using Interpolations
# using Optim
# using RCall
# R"library(stabledist)"

function checklevy()
    xr=-6.0:0.05:5.0
    mu=2.0
    sigma=0.5
    a=0.5 # Levy
    b=1.0 # Levy
    for pm in [0 1]
        if pm == 0
            c=sigma
            d=mu + sigma
        else
            c=sigma
            d=mu
        end
        dl=Levy(mu,sigma)
        yl=pdf.(dl,xr)
        ds=Stablep(a,b,c,d;p=pm)
        ys=pdf.(ds,xr)
        p=plot(title="Levy vs Stable, pm=$pm",legend=:outerbottom)
        p=plot!(p,xr,yl,label="Levy($mu,$sigma)")
        p=plot!(p,xr,ys,label="Sp($a,$b,$c,$d;p=$pm)")
        display(p)
        sleep(5)
    end
end

function checkcauchy()
    xr=-6.0:0.05:5.0
    mu=2.0
    sigma=0.5
    a=1.01 # Cauchy
    b=0.0
    c=sigma
    d=mu
    for pm in 2 # [0 1 2]
        dc=Cauchy(mu,sigma)
        yc=pdf.(dc,xr)
        ds=Stablep(a,b,c,d;p=pm)
        ys=pdf.(ds,xr)
        p=plot(title="Cauchy vs Stable, pm=$pm",legend=:outerbottom)
        p=plot!(p,xr,yc,label="Cauchy($mu,$sigma)")
        p=plot!(p,xr,ys,label="Sp($a,$b,$c,$d;p=$pm)")
        display(p)
        sleep(5)
    end
end

function checknormal()
    xr=-6.0:0.05:5.0
    mu=2.0
    sigma=0.5
    a=2.0 # Normal
    b=0.0
    c=sigma
    d=mu
    for pm in 2 # [0 1 2]
        dn=Normal(mu,sigma)
        yn=pdf.(dn,xr)
        ds=Stablep(a,b,c,d;p=pm)
        ys=pdf.(ds,xr)
        p=plot(title="Normal vs Stable, pm=$pm",legend=:outerbottom)
        p=plot!(p,xr,yn,label="N($mu,$sigma)")
        p=plot!(p,xr,ys,label="Sp($a,$b,$c,$d;p=$pm)")
        display(p)
        sleep(5)
    end
end


checklevy()
checkcauchy()
checknormal()

nothing
