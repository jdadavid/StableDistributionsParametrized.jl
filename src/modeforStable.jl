# """
#    mode(d::Stable) : compute mode for for Stable distribution
#    see Nolan's chapter 3.2, p.77 for _mm function
# """
function mode(d::Stable)
    α,β,γ,δ = params(d)
    d=0.0
    if !(α ≈ 1.0)
        d=δ+β*γ*tanpi(α/2)
    else
        d=δ+β*2/π*γ*log(γ)
    end
    return γ*_mm(α,β)+d
end
