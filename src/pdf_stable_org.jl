using SpecialFunctions
using QuadGK: quadgk

# integral representation from Nolan ch. 3
function pdf(d::Stable{T}, x::Real) where T
    α, β, σ, μ =  params(d)

    α == 2one(T) && return pdf(Normal(μ, √2σ),x)
    α == one(T) && β == zero(T) && return pdf(Cauchy(μ,σ),x)
    α == one(T)/2 && β == one(T) && return pdf(Levy(μ, σ), x)
    α == one(T)/2 && β == -one(T) && return pdf(Levy(-μ, σ), -x)

    w(v,c) = v*c > 36. ? 0.0 : v*exp(-c*v) # numerical truncation

    if α == one(T) 
        V₁(θ) = 2/π*(π/2+β*θ)/cos(θ) * exp((π/2+β*θ)*tan(θ)/β)

        x = (x-μ)/σ - 2/π*β*log(σ) # normalize to S(1,β,1,0)
        x < 0 && ( (x, β, μ) = (-x, -β, -μ) ) # reflection property

        I, _err = quadgk(θ -> w(V₁(θ),exp(-π*x/2β)), -π/2, π/2 ) 

        return 1/(2abs(β)*σ) * exp(-π*x/2β) * I
    else 
        V(θ) =(cos(α*θ₀))^(1/(α-1)) * (cos(θ)/sin(α*(θ₀+θ)))^(α/(α-1)) * cos(α*θ₀ + (α-1)*θ)/cos(θ)

        x = (x-μ)/σ # normalize to S(α,β,1,0)
        x < 0 && ( (x, β, μ) = (-x, -β, -μ) ) # reflection property

        θ₀ =  atan(β*tan(α*π/2))/α
        x ≈ 0. && return gamma(1+1/α)*cos(θ₀)*(cos(α*θ₀))^(1/α) / π
        
        I, _err =  quadgk(θ -> w(V(θ), x^(α/(α-1)) ), -θ₀, π/2)

        return α/σ * x^(1/(α-1)) / (π*abs(α-1)) * I
    end
end

