using Distributions
using StableDistributions
using Plots

include("pdf_stable.jl")
include("r_dpqrstable.jl")

function pb0re1()
    d=Stable(0.75,0.5)
    dpdf(x)=pdf(d,x)
    p=plot(title="pb0mwe : pb0re1 : Stable(0.75,0.5)")
    xr=-0.010:0.00001:0.02
    yd=dpdf.(xr)
    p=plot!(p,xr,yd,label="pdf")
    display(p)
end
function pb0mwe()
    d=Stable(0.75,0.5)
    x=0.005
    y=pdf(d,x)
    ry=r_dstable(x,0.75,0.5,    1.0, 0.0,  1)
    println("pb0mwe : Stable(0.75,0.5) : x=",x," pdf=",y," pdf_from_R=",ry)
end

@info "pb0mwe"; pb0mwe()
# @info "pb0re1"; pb0re1()
