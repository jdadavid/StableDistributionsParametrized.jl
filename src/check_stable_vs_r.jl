using Distributions
using StableDistributions
using StatsBase
using StatsPlots

include("r_dpqrstable.jl")

function checkvsrd()
    pm=1
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

    for c in gammar
        for d in deltar
            mmdy=0.0
            for a in alphar
                for b in betar
                    if a ≈ 2.0; b=0.0; end
                    # @show a,b,c,d; # sleep(2)
                    # @info "a=$a b=$b c=$c d=$d"
                    yr=r_dstabler.(xr,a,b,c,d,pm)
                    ds=Stable(a,b,c,d)
                    ys=pdf.(ds,xr)
                    dy=abs.(ys-yr)
                    mdy=maximum(dy)
                    # @show a,b,c,d,mdy
                    if mdy < mdyhuge; mmdy=max(mdy,mmdy); end
                    if mdy >= mdyhuge; @info "mdy huge $mdy for a=$a b=$b c=$c d=$d"; end # sleep(10);end
                    if doplot
                        p=plot(title="pdf",legend=:outerbottom)
                        p=plot!(p,xr,yr,label="R pdf a=$a b=$b c=$c d=$d")
                        p=plot!(p,xr,ys,label="J pdf a=$a b=$b c=$c d=$d")
                        display(p)
                        sleep(2)
                    end
                end # b
            end # a
            @show c,d,mmdy
        end # d
    end # c
end

function checkvsrp()
    pm=1
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
    for c in gammar
        for d in deltar
            mmdy=0.0
            for a in alphar
                for b in betar
                    if a ≈ 2.0; b=0.0; end
                    # @show a,b,c,d; # sleep(2)
                    # @info "a=$a b=$b c=$c d=$d"
                    yr=r_pstabler.(xr,a,b,c,d,pm)
                    ds=Stable(a,b,c,d)
                    ys=cdf.(ds,xr)
                    dy=abs.(ys-yr)
                    mdy=maximum(dy)
                    # @show a,b,c,d,mdy
                    if mdy < mdyhuge; mmdy=max(mdy,mmdy); end
                    if mdy >= mdyhuge; @info "mdy huge $mdy for a=$a b=$b c=$c d=$d"; end # sleep(10);end
                    if doplot
                        p=plot(title="cdf",legend=:outerbottom)
                        p=plot!(p,xr,yr,label="R cdf a=$a b=$b c=$c d=$d")
                        p=plot!(p,xr,ys,label="J cdf a=$a b=$b c=$c d=$d")
                        display(p)
                        sleep(2)
                    end
                end # b
            end # a
            @show c,d,mmdy
        end # d
    end # c
end

checkvsrd()





nothing
