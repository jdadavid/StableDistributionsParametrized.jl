using Distributions
using StableDistributions
if !isdefined(Main,:StableDistributionsp); include("StableDistributionsp.jl"); using ..StableDistributionsp; end
using StatsBase
using StatsPlots
# using Interpolations
# using Optim

include("r_dpqrstable.jl")

function checkvsrd()
alphaeps=0.001
alphar=vcat(alphaeps,collect(0.0:0.05:0.95),1.0-alphaeps,1.0,1.0+alphaeps,collect(1.05:0.05:2.0))
betar=0.:0.05:1.0
gammar=1.0
deltar=0.0
gammar=2.0
deltar=0.5
xr=-5.0:0.05:5.0

alphar=[0.3,0.5,1.0,1.3,2.0]
betar=[0.0,0.5,1.0]
# alphar=[1.3]
# betar=[0.5]
# alphar=[0.3]

mdyhuge=0.001
doplot=true
for pm in [0 1 2]
    for c in gammar
        for d in deltar
            mmdy=0.0
            for a in alphar
                for b in betar
                    if a ≈ 2.0; b=0.0; end
                    # @show pm,a,b,c,d; # sleep(2)
                    # @info "pm=$pm a=$a b=$b c=$c d=$d"
                    yr=r_dstabler.(xr,a,b,c,d,pm)
                    ds=Stablep(a,b,c,d;p=pm)
                    ys=pdf.(ds,xr)
                    dy=abs.(ys-yr)
                    mdy=maximum(dy)
                    # @show pm,a,b,c,d,mdy
                    if mdy < mdyhuge; mmdy=max(mdy,mmdy); end
                    if mdy >= mdyhuge; @info "mdy huge $mdy for pm=$pm a=$a b=$b c=$c d=$d"; end # sleep(10);end
                    if doplot
                        p=plot(title="pdf pm=$pm",legend=:outerbottom)
                        p=plot!(p,xr,yr,label="R pdf a=$a b=$b c=$c d=$d pm=$pm")
                        p=plot!(p,xr,ys,label="J pdf a=$a b=$b c=$c d=$d pm=$pm")
                        display(p)
                        sleep(2)
                    end
                end # b
            end # a
            @show pm,c,d,mmdy
        end # d
    end # c
end # pm
end

function checkvsrp()
alphaeps=0.001
alphar=vcat(alphaeps,collect(0.0:0.05:0.95),1.0-alphaeps,1.0,1.0+alphaeps,collect(1.05:0.05:2.0))
betar=0.:0.05:1.0
gammar=1.0
deltar=0.0
xr=-5.0:0.05:5.0

alphar=[0.3,0.5,1.0,1.3,2.0]
betar=[0.0,0.5,1.0]
# alphar=[1.3]
# betar=[0.5]
# alphar=[0.3]

mdyhuge=0.001
doplot=true
for pm in [0 1 2]
    for c in gammar
        for d in deltar
            mmdy=0.0
            for a in alphar
                for b in betar
                    if a ≈ 2.0; b=0.0; end
                    # @show pm,a,b,c,d; # sleep(2)
                    # @info "pm=$pm a=$a b=$b c=$c d=$d"
                    yr=r_pstabler.(xr,a,b,c,d,pm)
                    ds=Stablep(a,b,c,d;p=pm)
                    ys=cdf.(ds,xr)
                    dy=abs.(ys-yr)
                    mdy=maximum(dy)
                    # @show pm,a,b,c,d,mdy
                    if mdy < mdyhuge; mmdy=max(mdy,mmdy); end
                    if mdy >= mdyhuge; @info "mdy huge $mdy for pm=$pm a=$a b=$b c=$c d=$d"; end # sleep(10);end
                    if doplot
                        p=plot(title="cdf pm=$pm",legend=:outerbottom)
                        p=plot!(p,xr,yr,label="R cdf a=$a b=$b c=$c d=$d pm=$pm")
                        p=plot!(p,xr,ys,label="J cdf a=$a b=$b c=$c d=$d pm=$pm")
                        display(p)
                        sleep(2)
                    end
                end # b
            end # a
            @show pm,c,d,mmdy
        end # d
    end # c
end # pm
end

checkvsrd()





nothing
