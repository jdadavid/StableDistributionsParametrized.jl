using StableDistributions
using Distributions
using Interpolations
using Optim
using StatsBase
using StatsPlots

global _stablexdbg=false
global _stablexdbg1=false

"""
    stabledists_paramconv(pfrom,pto,α,β,γ,δ) : convert parametrization "pfrom" to parametrization "pto"
      see table 3.1 of Chapter 3.5 (p.97) of Univariate Stable Distributions - John P. Nolan (2020)
"""
function stabledists_paramconv(pfrom,pto,α,β,γ,δ)
    if pfrom == pto ; return α,β,γ,δ; end
    if !(pfrom in [0,1,2])
        error("stabledistributionsx : stabledists_paramconv : pfrom=$pfrom not in [0,1,2]")
    end
    if !(pto in [0,1,2])
        error("stabledistributionsx : stabledists_paramconv : pto=$pto not in [0,1,2]")
    end
    (a,b,c,d)=(α,β,γ,δ)

    if     pfrom == 0
        if     pto == 1
            if !(α ≈ 1.0)
                d=δ-β*γ*tanpi(α/2)
            else
                d=δ-β*2/π*γ*log(γ)
            end
        elseif pto == 2
            c=γ/(α^(-1.0/α))
            d=δ+γ*mm(α,β)
        end
    elseif pfrom == 1
        if     pto == 0
            if !(α ≈ 1.0)
                d=δ+β*γ*tanpi(α/2)
            else
                d=δ+β*2/π*γ*log(γ)
            end
        elseif pto == 2
            c=γ/(α^(-1.0/α))
            if !(α ≈ 1.0)
                d=δ+γ*(β*tanpi(α/2)+mm(α,β))
            else
                d=δ+c*(β*2/π*log(c)+mm(α,β))
            end
        end
    elseif pfrom == 2
        c=α^(-1.0/α)*γ
        if     pto == 0
            d=δ-c*mm(α,β)
        elseif pto == 1
            if !(α ≈ 1.0)
                d=δ-c*(β*tanpi(α/2)+mm(α,β))
            else
                d=δ-γ*(β*2/π*log(γ)+mm(α,β))
            end
        end
    end

    _stablexdbg && @info "pfrom=($α,$β,$γ,$δ)) => pto=($a,$b,$c,$d))"
    _stablexdbg1 && @info "pfrom=($α,$β,$γ,$δ)) => pto=($a,$b,$c,$d))"
    return (a,b,c,d)
end

"""
    mmv1(α,β) : compute mode for Stable distribution(α,β) with parametrization 0
      see Appendix C of Univariate Stable Distributions - John P. Nolan (2020) for data
      function is using internal (JD's) linear interpolation between Nolan's tabulated data
"""
function mmv1(α,β)
    if !(0.01 <= α <= 2.0)
        error("stabledistributionsx : mmv1 : α=$α not inside ]0.01,2.0]")
    end
    if !(-1.0 <= β <= 1.0)
        error("stabledistributionsx : mmv1 : β=$β not inside [-1.0,1.0]")
    end
    alpha=α
    beta=β
    if β<0 ; betaneg=true; beta=-β; else; betaneg=false; end
    
# from Appendix C of Univariate Stable Distributions - John P. Nolan (2020)
# Line for α=0 contains invalid values, but is used for computing index for interpolation
 #    β =   0.0      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0   #  α
    A  = [0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000 # 0.00 INVALID?
          0.00000 -0.00787 -0.01574 -0.02361 -0.03148 -0.03935 -0.04722 -0.05509 -0.06296 -0.07083 -0.07870 # 0.05
          0.00000 -0.01584 -0.03168 -0.04752 -0.06335 -0.07919 -0.09503 -0.11087 -0.12671 -0.14255 -0.15338 # 0.10
          0.00000 -0.02401 -0.04802 -0.07202 -0.09603 -0.12004 -0.14405 -0.16806 -0.19206 -0.21607 -0.23508 # 0.15
          0.00000 -0.03249 -0.06498 -0.09748 -0.12997 -0.16246 -0.19495 -0.22744 -0.25994 -0.28743 -0.31992 # 0.20
          0.00000 -0.04142 -0.08284 -0.12426 -0.16569 -0.20711 -0.24353 -0.28495 -0.32637 -0.36779 -0.40921 # 0.25
          0.00000 -0.05095 -0.10191 -0.14786 -0.19881 -0.24976 -0.29986 -0.34879 -0.39731 -0.44539 -0.49299 # 0.30
          0.00000 -0.06128 -0.11756 -0.17705 -0.23479 -0.29177 -0.34792 -0.40319 -0.45751 -0.51082 -0.56307 # 0.35
          0.00000 -0.06765 -0.13452 -0.20004 -0.26418 -0.32686 -0.38800 -0.44754 -0.50540 -0.56153 -0.61587 # 0.40
          0.00000 -0.07454 -0.14754 -0.21852 -0.28730 -0.35380 -0.41794 -0.47967 -0.53894 -0.59572 -0.64995 # 0.45
          0.00000 -0.07955 -0.15704 -0.23174 -0.30340 -0.37193 -0.43729 -0.49945 -0.55841 -0.61415 -0.66667 # 0.50
          0.00000 -0.08267 -0.16288 -0.23963 -0.31258 -0.38163 -0.44678 -0.50804 -0.56544 -0.61902 -0.66880 # 0.55
          0.00000 -0.08399 -0.16527 -0.24259 -0.31549 -0.38388 -0.44776 -0.50719 -0.56223 -0.61297 -0.65948 # 0.60
          0.00000 -0.08372 -0.16463 -0.24127 -0.31307 -0.37990 -0.44177 -0.49876 -0.55099 -0.59857 -0.64162 # 0.65
          0.00000 -0.08213 -0.16147 -0.23642 -0.30630 -0.37092 -0.43030 -0.48452 -0.53372 -0.57806 -0.61768 # 0.70
          0.00000 -0.07949 -0.15631 -0.22878 -0.29613 -0.35810 -0.41468 -0.46597 -0.51211 -0.55327 -0.58963 # 0.75
          0.00000 -0.07606 -0.14964 -0.21904 -0.28341 -0.34243 -0.39605 -0.44436 -0.48750 -0.52565 -0.55899 # 0.80
          0.00000 -0.07208 -0.14188 -0.20777 -0.26886 -0.32474 -0.37534 -0.42070 -0.46097 -0.49631 -0.52692 # 0.85
          0.00000 -0.06772 -0.13339 -0.19549 -0.25308 -0.30572 -0.35328 -0.39577 -0.43332 -0.46609 -0.49424 # 0.90
          0.00000 -0.06314 -0.12448 -0.18259 -0.23656 -0.28592 -0.33047 -0.37020 -0.40519 -0.43559 -0.46156 # 0.95
          0.00000 -0.05847 -0.11537 -0.16940 -0.21970 -0.26577 -0.30736 -0.34443 -0.37704 -0.40528 -0.42931 # 1.00
          0.00000 -0.05380 -0.10625 -0.15619 -0.20281 -0.24561 -0.28431 -0.31884 -0.34920 -0.37547 -0.39778 # 1.05
          0.00000 -0.04921 -0.09725 -0.14315 -0.18613 -0.22571 -0.26160 -0.29367 -0.32192 -0.34639 -0.36717 # 1.10
          0.00000 -0.04474 -0.08850 -0.13042 -0.16983 -0.20626 -0.23941 -0.26914 -0.29539 -0.31818 -0.33758 # 1.15
          0.00000 -0.04043 -0.08005 -0.11812 -0.15406 -0.18742 -0.21791 -0.24537 -0.26972 -0.29094 -0.30909 # 1.20
          0.00000 -0.03632 -0.07197 -0.10633 -0.13891 -0.16929 -0.19720 -0.22247 -0.24499 -0.26474 -0.28173 # 1.25
          0.00000 -0.03242 -0.06429 -0.09511 -0.12444 -0.15195 -0.17736 -0.20051 -0.22127 -0.23961 -0.25550 # 1.30
          0.00000 -0.02873 -0.05703 -0.08448 -0.11072 -0.13545 -0.15845 -0.17953 -0.19859 -0.21555 -0.23040 # 1.35
          0.00000 -0.02527 -0.05020 -0.07446 -0.09775 -0.11983 -0.14049 -0.15957 -0.17696 -0.19258 -0.20639 # 1.40
          0.00000 -0.02204 -0.04381 -0.06507 -0.08556 -0.10510 -0.12350 -0.14063 -0.15639 -0.17068 -0.18347 # 1.45
          0.00000 -0.01903 -0.03786 -0.05629 -0.07415 -0.09126 -0.10750 -0.12273 -0.13687 -0.14984 -0.16159 # 1.50
          0.00000 -0.01624 -0.03233 -0.04813 -0.06350 -0.07832 -0.09248 -0.10587 -0.11842 -0.13006 -0.14073 # 1.55
          0.00000 -0.01366 -0.02722 -0.04058 -0.05362 -0.06627 -0.07843 -0.09003 -0.10101 -0.11131 -0.12088 # 1.60
          0.00000 -0.01130 -0.02252 -0.03361 -0.04448 -0.05508 -0.06534 -0.07522 -0.08465 -0.09360 -0.10202 # 1.65
          0.00000 -0.00913 -0.01822 -0.02721 -0.03607 -0.04475 -0.05321 -0.06142 -0.06933 -0.07692 -0.08415 # 1.70
          0.00000 -0.00716 -0.01429 -0.02137 -0.02837 -0.03526 -0.04202 -0.04863 -0.05505 -0.06128 -0.06729 # 1.75
          0.00000 -0.00538 -0.01074 -0.01607 -0.02136 -0.02660 -0.03176 -0.03684 -0.04183 -0.04670 -0.05146 # 1.80
          0.00000 -0.00377 -0.00754 -0.01130 -0.01503 -0.01874 -0.02243 -0.02607 -0.02967 -0.03323 -0.03673 # 1.85
          0.00000 -0.00235 -0.00469 -0.00703 -0.00937 -0.01170 -0.01401 -0.01632 -0.01861 -0.02089 -0.02315 # 1.90
          0.00000 -0.00109 -0.00218 -0.00327 -0.00436 -0.00545 -0.00653 -0.00762 -0.00870 -0.00978 -0.01086 # 1.95
          0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000]# 2.00
    _stablexdbg && @show size(A)

    valmm=0.0

    ialpha=round(Int,20*alpha,RoundDown)
    alphar=ialpha/20.0
    alphaf=alpha-alphar
    falpha=20*alphaf
    ialpha=1+ialpha
    _stablexdbg && @show alpha,alphar,alphaf,ialpha,falpha
    if (ialpha<1) || (ialpha>size(A,1)) ; error("stabledistributionsx : mmo1de : ialpha=$ialpha not in range [1,$(size(A,1))]"); end
    if (falpha<0.0) || (falpha>1.0)     ; error("stabledistributionsx : mmo1de : falpha=$falpha not in range [0.0,1.0]"); end

    ibeta=round(Int,10*beta,RoundDown)
    betar=ibeta/10.0
    betaf=beta-betar
    fbeta=10*betaf
    ibeta=1+ibeta
    _stablexdbg && @show beta,betar,betaf,ibeta,fbeta
    if (ibeta<1)   || (ibeta>size(A,2)) ; error("stabledistributionsx : mmo1de : ibeta=$ibeta not in range [1,$(size(A,2))]"); end
    if (fbeta<0.0) || (fbeta>1.0)       ; error("stabledistributionsx : mmo1de : fbeta=$fbeta not in range [0.0,1.0]"); end

    alphamax = (ialpha == size(A,1))
    betamax  = (ibeta  == size(A,2))
    if alphamax && betamax
        valmm=A[ialpha,ibeta]
        _stablexdbg && @show "corner",valmm
    elseif alphamax
        valb0=A[ialpha,ibeta]
        valb1=A[ialpha,ibeta+1]
        valmm=valb0+fbeta*(valb1-valb0)
        _stablexdbg && @show "alphamax",valb0,valb1,fbeta,valmm
    elseif betamax
        vala0=A[ialpha  ,ibeta]
        vala1=A[ialpha+1,ibeta]
        valmm=vala0+falpha*(vala1-vala0)
        _stablexdbg && @show "betamax",vala0,vala1,falpha,valmm
    else
        vala0=A[ialpha  ,ibeta] + fbeta*(A[ialpha  ,ibeta+1]-A[ialpha  ,ibeta])
        vala1=A[ialpha+1,ibeta] + fbeta*(A[ialpha+1,ibeta+1]-A[ialpha+1,ibeta])
        valmm=vala0+falpha*(vala1-vala0)
        _stablexdbg && @show "std",vala0,vala1,falpha,fbeta,valmm
    end

    if betaneg; valmm=-valmm; end
    return valmm
end

"""
    mmv2(α,β) : compute mode for Stable distribution(α,β) with parametrization 0
      see Appendix C of Univariate Stable Distributions - John P. Nolan (2020) for data
      function is using Julia's interpolation.jl between Nolan's tabulated data
"""
function mmv2(α,β;cubic=true)
    if !(0.01 <= α <= 2.0)
        error("stabledistributionsx : mmv2 : α=$α not inside ]0.01,2.0]")
    end
    if !(-1.0 <= β <= 1.0)
        error("stabledistributionsx : mmv2 : β=$β not inside [-1.0,1.0]")
    end
    alpha=α
    beta=β
    if β<0 ; betaneg=true; beta=-β; else; betaneg=false; end

    xs=0.0:0.05:2.0 # alpha range
    ys=0.0:0.1:1.0  # beta  range
# from Appendix C of Univariate Stable Distributions - John P. Nolan (2020)
# Line for α=0 contains invalid values, but is used for computing index for interpolation
 #    β =   0.0      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0   #  α
    A  = [0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000 # 0.00 INVALID?
          0.00000 -0.00787 -0.01574 -0.02361 -0.03148 -0.03935 -0.04722 -0.05509 -0.06296 -0.07083 -0.07870 # 0.05
          0.00000 -0.01584 -0.03168 -0.04752 -0.06335 -0.07919 -0.09503 -0.11087 -0.12671 -0.14255 -0.15338 # 0.10
          0.00000 -0.02401 -0.04802 -0.07202 -0.09603 -0.12004 -0.14405 -0.16806 -0.19206 -0.21607 -0.23508 # 0.15
          0.00000 -0.03249 -0.06498 -0.09748 -0.12997 -0.16246 -0.19495 -0.22744 -0.25994 -0.28743 -0.31992 # 0.20
          0.00000 -0.04142 -0.08284 -0.12426 -0.16569 -0.20711 -0.24353 -0.28495 -0.32637 -0.36779 -0.40921 # 0.25
          0.00000 -0.05095 -0.10191 -0.14786 -0.19881 -0.24976 -0.29986 -0.34879 -0.39731 -0.44539 -0.49299 # 0.30
          0.00000 -0.06128 -0.11756 -0.17705 -0.23479 -0.29177 -0.34792 -0.40319 -0.45751 -0.51082 -0.56307 # 0.35
          0.00000 -0.06765 -0.13452 -0.20004 -0.26418 -0.32686 -0.38800 -0.44754 -0.50540 -0.56153 -0.61587 # 0.40
          0.00000 -0.07454 -0.14754 -0.21852 -0.28730 -0.35380 -0.41794 -0.47967 -0.53894 -0.59572 -0.64995 # 0.45
          0.00000 -0.07955 -0.15704 -0.23174 -0.30340 -0.37193 -0.43729 -0.49945 -0.55841 -0.61415 -0.66667 # 0.50
          0.00000 -0.08267 -0.16288 -0.23963 -0.31258 -0.38163 -0.44678 -0.50804 -0.56544 -0.61902 -0.66880 # 0.55
          0.00000 -0.08399 -0.16527 -0.24259 -0.31549 -0.38388 -0.44776 -0.50719 -0.56223 -0.61297 -0.65948 # 0.60
          0.00000 -0.08372 -0.16463 -0.24127 -0.31307 -0.37990 -0.44177 -0.49876 -0.55099 -0.59857 -0.64162 # 0.65
          0.00000 -0.08213 -0.16147 -0.23642 -0.30630 -0.37092 -0.43030 -0.48452 -0.53372 -0.57806 -0.61768 # 0.70
          0.00000 -0.07949 -0.15631 -0.22878 -0.29613 -0.35810 -0.41468 -0.46597 -0.51211 -0.55327 -0.58963 # 0.75
          0.00000 -0.07606 -0.14964 -0.21904 -0.28341 -0.34243 -0.39605 -0.44436 -0.48750 -0.52565 -0.55899 # 0.80
          0.00000 -0.07208 -0.14188 -0.20777 -0.26886 -0.32474 -0.37534 -0.42070 -0.46097 -0.49631 -0.52692 # 0.85
          0.00000 -0.06772 -0.13339 -0.19549 -0.25308 -0.30572 -0.35328 -0.39577 -0.43332 -0.46609 -0.49424 # 0.90
          0.00000 -0.06314 -0.12448 -0.18259 -0.23656 -0.28592 -0.33047 -0.37020 -0.40519 -0.43559 -0.46156 # 0.95
          0.00000 -0.05847 -0.11537 -0.16940 -0.21970 -0.26577 -0.30736 -0.34443 -0.37704 -0.40528 -0.42931 # 1.00
          0.00000 -0.05380 -0.10625 -0.15619 -0.20281 -0.24561 -0.28431 -0.31884 -0.34920 -0.37547 -0.39778 # 1.05
          0.00000 -0.04921 -0.09725 -0.14315 -0.18613 -0.22571 -0.26160 -0.29367 -0.32192 -0.34639 -0.36717 # 1.10
          0.00000 -0.04474 -0.08850 -0.13042 -0.16983 -0.20626 -0.23941 -0.26914 -0.29539 -0.31818 -0.33758 # 1.15
          0.00000 -0.04043 -0.08005 -0.11812 -0.15406 -0.18742 -0.21791 -0.24537 -0.26972 -0.29094 -0.30909 # 1.20
          0.00000 -0.03632 -0.07197 -0.10633 -0.13891 -0.16929 -0.19720 -0.22247 -0.24499 -0.26474 -0.28173 # 1.25
          0.00000 -0.03242 -0.06429 -0.09511 -0.12444 -0.15195 -0.17736 -0.20051 -0.22127 -0.23961 -0.25550 # 1.30
          0.00000 -0.02873 -0.05703 -0.08448 -0.11072 -0.13545 -0.15845 -0.17953 -0.19859 -0.21555 -0.23040 # 1.35
          0.00000 -0.02527 -0.05020 -0.07446 -0.09775 -0.11983 -0.14049 -0.15957 -0.17696 -0.19258 -0.20639 # 1.40
          0.00000 -0.02204 -0.04381 -0.06507 -0.08556 -0.10510 -0.12350 -0.14063 -0.15639 -0.17068 -0.18347 # 1.45
          0.00000 -0.01903 -0.03786 -0.05629 -0.07415 -0.09126 -0.10750 -0.12273 -0.13687 -0.14984 -0.16159 # 1.50
          0.00000 -0.01624 -0.03233 -0.04813 -0.06350 -0.07832 -0.09248 -0.10587 -0.11842 -0.13006 -0.14073 # 1.55
          0.00000 -0.01366 -0.02722 -0.04058 -0.05362 -0.06627 -0.07843 -0.09003 -0.10101 -0.11131 -0.12088 # 1.60
          0.00000 -0.01130 -0.02252 -0.03361 -0.04448 -0.05508 -0.06534 -0.07522 -0.08465 -0.09360 -0.10202 # 1.65
          0.00000 -0.00913 -0.01822 -0.02721 -0.03607 -0.04475 -0.05321 -0.06142 -0.06933 -0.07692 -0.08415 # 1.70
          0.00000 -0.00716 -0.01429 -0.02137 -0.02837 -0.03526 -0.04202 -0.04863 -0.05505 -0.06128 -0.06729 # 1.75
          0.00000 -0.00538 -0.01074 -0.01607 -0.02136 -0.02660 -0.03176 -0.03684 -0.04183 -0.04670 -0.05146 # 1.80
          0.00000 -0.00377 -0.00754 -0.01130 -0.01503 -0.01874 -0.02243 -0.02607 -0.02967 -0.03323 -0.03673 # 1.85
          0.00000 -0.00235 -0.00469 -0.00703 -0.00937 -0.01170 -0.01401 -0.01632 -0.01861 -0.02089 -0.02315 # 1.90
          0.00000 -0.00109 -0.00218 -0.00327 -0.00436 -0.00545 -0.00653 -0.00762 -0.00870 -0.00978 -0.01086 # 1.95
          0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000]# 2.00

    valmm=0.0

    if !cubic
        # linear interpolation
        interp_linear = linear_interpolation((xs, ys), A)
        valmm=interp_linear(alpha,beta)
    else
        # cubic spline interpolation
        interp_cubic = cubic_spline_interpolation((xs, ys), A)
        valmm=interp_cubic(alpha,beta)
    end

    if betaneg; valmm=-valmm; end
    return valmm
end


"""
    mmv3(α,β) : compute mode for Stable distribution(α,β) with parametrization 0
      see Appendix C of Univariate Stable Distributions - John P. Nolan (2020) for data
      function is using  linear interpolation (from StableDistribution src) between Nolan's tabulated data
"""
function mmv3(α,β)
    if !(0.01 <= α <= 2.0)
        error("stabledistributionsx : mmv1 : α=$α not inside ]0.01,2.0]")
    end
    if !(-1.0 <= β <= 1.0)
        error("stabledistributionsx : mmv1 : β=$β not inside [-1.0,1.0]")
    end
    alpha=α
    beta=β
    if β<0 ; betaneg=true; beta=-β; else; betaneg=false; end

    xs=0.0:0.05:2.0 # alpha range
    ys=0.0:0.1:1.0  # beta  range

# from Appendix C of Univariate Stable Distributions - John P. Nolan (2020)
# Line for α=0 contains invalid values, but is used for computing index for interpolation
 #    β =   0.0      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0   #  α
    A  = [0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000 # 0.00 INVALID?
          0.00000 -0.00787 -0.01574 -0.02361 -0.03148 -0.03935 -0.04722 -0.05509 -0.06296 -0.07083 -0.07870 # 0.05
          0.00000 -0.01584 -0.03168 -0.04752 -0.06335 -0.07919 -0.09503 -0.11087 -0.12671 -0.14255 -0.15338 # 0.10
          0.00000 -0.02401 -0.04802 -0.07202 -0.09603 -0.12004 -0.14405 -0.16806 -0.19206 -0.21607 -0.23508 # 0.15
          0.00000 -0.03249 -0.06498 -0.09748 -0.12997 -0.16246 -0.19495 -0.22744 -0.25994 -0.28743 -0.31992 # 0.20
          0.00000 -0.04142 -0.08284 -0.12426 -0.16569 -0.20711 -0.24353 -0.28495 -0.32637 -0.36779 -0.40921 # 0.25
          0.00000 -0.05095 -0.10191 -0.14786 -0.19881 -0.24976 -0.29986 -0.34879 -0.39731 -0.44539 -0.49299 # 0.30
          0.00000 -0.06128 -0.11756 -0.17705 -0.23479 -0.29177 -0.34792 -0.40319 -0.45751 -0.51082 -0.56307 # 0.35
          0.00000 -0.06765 -0.13452 -0.20004 -0.26418 -0.32686 -0.38800 -0.44754 -0.50540 -0.56153 -0.61587 # 0.40
          0.00000 -0.07454 -0.14754 -0.21852 -0.28730 -0.35380 -0.41794 -0.47967 -0.53894 -0.59572 -0.64995 # 0.45
          0.00000 -0.07955 -0.15704 -0.23174 -0.30340 -0.37193 -0.43729 -0.49945 -0.55841 -0.61415 -0.66667 # 0.50
          0.00000 -0.08267 -0.16288 -0.23963 -0.31258 -0.38163 -0.44678 -0.50804 -0.56544 -0.61902 -0.66880 # 0.55
          0.00000 -0.08399 -0.16527 -0.24259 -0.31549 -0.38388 -0.44776 -0.50719 -0.56223 -0.61297 -0.65948 # 0.60
          0.00000 -0.08372 -0.16463 -0.24127 -0.31307 -0.37990 -0.44177 -0.49876 -0.55099 -0.59857 -0.64162 # 0.65
          0.00000 -0.08213 -0.16147 -0.23642 -0.30630 -0.37092 -0.43030 -0.48452 -0.53372 -0.57806 -0.61768 # 0.70
          0.00000 -0.07949 -0.15631 -0.22878 -0.29613 -0.35810 -0.41468 -0.46597 -0.51211 -0.55327 -0.58963 # 0.75
          0.00000 -0.07606 -0.14964 -0.21904 -0.28341 -0.34243 -0.39605 -0.44436 -0.48750 -0.52565 -0.55899 # 0.80
          0.00000 -0.07208 -0.14188 -0.20777 -0.26886 -0.32474 -0.37534 -0.42070 -0.46097 -0.49631 -0.52692 # 0.85
          0.00000 -0.06772 -0.13339 -0.19549 -0.25308 -0.30572 -0.35328 -0.39577 -0.43332 -0.46609 -0.49424 # 0.90
          0.00000 -0.06314 -0.12448 -0.18259 -0.23656 -0.28592 -0.33047 -0.37020 -0.40519 -0.43559 -0.46156 # 0.95
          0.00000 -0.05847 -0.11537 -0.16940 -0.21970 -0.26577 -0.30736 -0.34443 -0.37704 -0.40528 -0.42931 # 1.00
          0.00000 -0.05380 -0.10625 -0.15619 -0.20281 -0.24561 -0.28431 -0.31884 -0.34920 -0.37547 -0.39778 # 1.05
          0.00000 -0.04921 -0.09725 -0.14315 -0.18613 -0.22571 -0.26160 -0.29367 -0.32192 -0.34639 -0.36717 # 1.10
          0.00000 -0.04474 -0.08850 -0.13042 -0.16983 -0.20626 -0.23941 -0.26914 -0.29539 -0.31818 -0.33758 # 1.15
          0.00000 -0.04043 -0.08005 -0.11812 -0.15406 -0.18742 -0.21791 -0.24537 -0.26972 -0.29094 -0.30909 # 1.20
          0.00000 -0.03632 -0.07197 -0.10633 -0.13891 -0.16929 -0.19720 -0.22247 -0.24499 -0.26474 -0.28173 # 1.25
          0.00000 -0.03242 -0.06429 -0.09511 -0.12444 -0.15195 -0.17736 -0.20051 -0.22127 -0.23961 -0.25550 # 1.30
          0.00000 -0.02873 -0.05703 -0.08448 -0.11072 -0.13545 -0.15845 -0.17953 -0.19859 -0.21555 -0.23040 # 1.35
          0.00000 -0.02527 -0.05020 -0.07446 -0.09775 -0.11983 -0.14049 -0.15957 -0.17696 -0.19258 -0.20639 # 1.40
          0.00000 -0.02204 -0.04381 -0.06507 -0.08556 -0.10510 -0.12350 -0.14063 -0.15639 -0.17068 -0.18347 # 1.45
          0.00000 -0.01903 -0.03786 -0.05629 -0.07415 -0.09126 -0.10750 -0.12273 -0.13687 -0.14984 -0.16159 # 1.50
          0.00000 -0.01624 -0.03233 -0.04813 -0.06350 -0.07832 -0.09248 -0.10587 -0.11842 -0.13006 -0.14073 # 1.55
          0.00000 -0.01366 -0.02722 -0.04058 -0.05362 -0.06627 -0.07843 -0.09003 -0.10101 -0.11131 -0.12088 # 1.60
          0.00000 -0.01130 -0.02252 -0.03361 -0.04448 -0.05508 -0.06534 -0.07522 -0.08465 -0.09360 -0.10202 # 1.65
          0.00000 -0.00913 -0.01822 -0.02721 -0.03607 -0.04475 -0.05321 -0.06142 -0.06933 -0.07692 -0.08415 # 1.70
          0.00000 -0.00716 -0.01429 -0.02137 -0.02837 -0.03526 -0.04202 -0.04863 -0.05505 -0.06128 -0.06729 # 1.75
          0.00000 -0.00538 -0.01074 -0.01607 -0.02136 -0.02660 -0.03176 -0.03684 -0.04183 -0.04670 -0.05146 # 1.80
          0.00000 -0.00377 -0.00754 -0.01130 -0.01503 -0.01874 -0.02243 -0.02607 -0.02967 -0.03323 -0.03673 # 1.85
          0.00000 -0.00235 -0.00469 -0.00703 -0.00937 -0.01170 -0.01401 -0.01632 -0.01861 -0.02089 -0.02315 # 1.90
          0.00000 -0.00109 -0.00218 -0.00327 -0.00436 -0.00545 -0.00653 -0.00762 -0.00870 -0.00978 -0.01086 # 1.95
          0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000]# 2.00

    valmm=0.0

    function interpolate2D(x, y, xs, ys, M)
        lx, ly = length(xs), length(ys)
        i, j  = findlast(xs .<= x), findlast(ys .<= y)
        i == lx && j < ly && return M[end, j]
        j == ly && i < lx && return M[i, end]
        i == lx && j == ly && return M[end,end]
        x₁, x₂, y₁, y₂ = xs[i], xs[i+1], ys[j], ys[j+1]
        m₁₁, m₁₂, m₂₁, m₂₂ = M[i, j], M[i, j+1], M[i+1, j], M[i+1, j+1]
        xmy = (x₂-x) * (m₁₁*(y₂-y) + m₁₂*(y-y₁)) + (x-x₁) * (m₂₁*(y₂-y) + m₂₂*(y-y₁))
        return xmy / ( (x₂-x₁)*(y₂-y₁) )
    end

    valmm=interpolate2D(alpha,beta,xs,ys,A)

    if betaneg; valmm=-valmm; end
    return valmm
end

"""
    mm(α,β) : compute mode for Stable distribution(α,β) with parametrization 0
      see Appendix C of Univariate Stable Distributions - John P. Nolan (2020) for data
"""
function mm(α,β)
    return mmv3(α,β)
end

"""
   mm0(α,β,γ,δ) : compute mode for for Stable distribution(α,β,γ,δ) with parametrization 0
   see Nolan's chapter 3.2, p.77
"""
function mm0(α,β,γ,δ)
    return γ*mm(α,β)+δ
end

"""
   mm0(α,β) : compute mode for for Stable distribution(α,β) with parametrization 0 (same as Nolan's mm)
"""
function mm0(α,β)
    return mm0(α,β,1.0,0.0)
end

"""
   mm0(α) : compute mode for for Stable distribution(α) with parametrization 0 (same as Nolan's mm)
"""
function mm0(α)
    return mm0(α,0.0,1.1,0.0)
end


"""
   mm1(α,β,γ,δ) : compute mode for for Stable distribution(α,β,γ,δ) with parametrization 1
   see Nolan's chapter 3.2, p.77
"""
function mm1(α,β,γ,δ)
    d=0.0
    if !(α ≈ 1.0)
        d=δ+β*γ*tanpi(α/2)
    else
        d=δ+β*2/π*γ*log(γ)
    end
    return γ*mm(α,β)+d
end

"""
   mm1(α,β) : compute mode for for Stable distribution(α,β) with parametrization 1
"""
function mm1(α,β)
    return mm1(α,β,1.0,0.0)
end

"""
   mm1(α) : compute mode for for Stable distribution(α) with parametrization 1
"""
function mm1(α)
    return mm1(α,0.0,1.1,0.0)
end


function mode(d::Stable)
    return mm1(params(d)...)
end

"""
   mm1old(α,β) : compute mode for for Stable distribution(α,β) with parametrization 1
"""
function mm1(α,β)
    (a,b,c,d)=stabledists_paramconv(1,0,α,β,1,0)
    _stablexdbg && @info "mm1 : shift by $d"
    _stablexdbg1 && @info "mm1 : shift by $d"
    return mm(a,b)+d # to be checked
end

"""
    mma1(α,β) : compute approx mode for Stable distribution(α,β) with parametrization 1
      see StableDistributions source
"""
function mma1(α,β)
    d=StableDistributions.Stable(α,β)
    return StableDistributions.appr_mode(d)
end

"""
    mma0(α,β) : compute approx mode for Stable distribution(α,β) with parametrization 0
      see StableDistributions source and see Nolan's chapter 3.2, p.77
"""
function mma0(α,β)
    γ=1.0; δ=0.0
    d=0.0
    if !(α ≈ 1.0)
        d=δ+β*γ*tanpi(α/2)
    else
        d=δ+β*2/π*γ*log(γ)
    end
    return mma1(α,β)-d
end

"""
    mmo1(α,β) : compute mode for Stable distribution(α,β) with parametrization 1
      use Optim.optimization
"""
function mmo1(α,β)
    atol=sqrt(eps(1.0))
    rtol=sqrt(eps(1.0))
    if isapprox(β,0.0; atol=atol,rtol=rtol); return 0.0;end # optimize fails if β==0 and α ≈ 1
    d=StableDistributions.Stable(α,β)
    g(x) = -pdf(d,x)
    res=Optim.optimize(g,-15.0,15.0)
    _stablexdbg && @show α,β,res
    return Optim.minimizer(res)
end

"""
    mmo1(α,β) : compute mode for Stable distribution(α,β) with parametrization 1
      use Optim.optimization
"""
function mmo0(α,β)
    γ=1.0; δ=0.0
    d=0.0
    if !(α ≈ 1.0)
        d=δ+β*γ*tanpi(α/2)
    else
        d=δ+β*2/π*γ*log(γ)
    end
    return mmo1(α,β)-d
end

function printvt()
    _stablexdbg=false
    ar=0.05:0.025:2.0
    for a in ar
        br=0.0:0.05:1.0
        # br=0.5
        println()
        if a==2.0; br=0.0;end
        for b in br
            print("a=$a b=$b : ")
            vt=mm(a,b)
            println("vt=$vt")
            #vo=mmo1(a,b)
            #dva=abs(vo-vt)
            #dvr=dva/max(vt,vo,1.e-10)
            #println("vt=$vt vo=$vo; dva=$dva dvr=$dvr")
        end
    end
end
# printvt()

function printdvar()
    _stablexdbg=false
    ar=0.1:0.1:2.0
    for a in ar
        println()
        br=0.0:0.1:1.0
        if a==2.0; br=0.0;end
        for b in br
            print("a=$a b=$b : ")
            vt=mm(a,b)
            vo=mmo1(a,b)
            dva=abs(vo-vt)
            dvr=dva/max(vt,vo,1.e-10)
            println("vt=$vt vo=$vo; dva=$dva dvr=$dvr")
        end
    end
end
# printdvar()

function plotmm()
    _stablexdbg=false
    p=plot(title="mm")
    ar=0.05:0.05:2.0
    for b in 0.0:0.2:1.0
        y=mm.(ar,b)
        p=plot!(ar,y,label="b=$b")
    end
    display(p)
end
function plotmma0()
    _stablexdbg=false
    ar=0.05:0.05:1.99
    b=1.0
    p2=plot(title="mm vs mma, b=$b")
    ymm0=mm0.(ar,b)
    ymma0=mma0.(ar,b)
    p2=plot!(ar,ymm0,label="mm0")
    p2=plot!(ar,ymma0,label="mma0")
    display(p2)
end
function plotmmao0()
    _stablexdbg=false
    ar=0.05:0.05:1.99
    b=1.0
    p2=plot(title="mm vs mma, b=1.0")
    ymm0=mm0.(ar,b)
    ymma0=mma0.(ar,b)
    ymmo0=mmo0.(ar,b)
    p2=plot!(ar,ymm0,label="mm0")
    p2=plot!(ar,ymma0,label="mma0")
    p2=plot!(ar,ymmo0,label="mmo0")
    display(p2)
end
function plotmmo0()
    _stablexdbg=false
    p=plot(title="mmo0")
    ar=0.05:0.05:1.99
    br=0.0:0.05:1.0
    # @show ar,br
    for b in br
        @show ar,b
        p=plot(title="mmo0")
        y=mm0.(ar,b)
        p=plot(ar,y,label="mm0  b=$b")
        display(p)
        # y=mmo0.(ar,b)
        y=similar(collect(ar))
        for (i,a) in pairs(ar)
            # @show a,b
            y[i]=mmo0(a,b)
            # @show i,y[i],a,b
        end
        p=plot!(ar,y,label="mmo0 b=$b")
        display(p)
        # sleep(0.5)
    end
    p=plot!()
    display(p)
end

function checkstabledists_paramconv()
    for (pfrom,pto) in ((0,1),(1,0),(0,2),(2,0),(1,2),(2,1))
        println(); @show pfrom,pto; println()
        ra=0.01:0.025:2.0
        rb=0.0:0.05:1.0
        c=1.0
        d=0.0
        epss=sqrt(eps(1.0))
        atolx=epss
        rtolx=1.0e-5
        atolc=atolx
        rtolc=rtolx
        atold=atolx
        rtold=rtolx
        for a in ra
            for b in rb
                (an,bn,cn,dn)=stabledists_paramconv(pfrom,pto,a,b,c,d)
                (aa,ba,ca,da)=stabledists_paramconv(pto,pfrom,an,bn,cn,dn)
                s=""
                if aa!=a;s="$s bad a=$a aa=$aa"; end
                if ba!=b;s="$s bad b=$b ba=$ba"; end
                if !isapprox(c,ca;atol=atolc,rtol=rtolc); s="$s bad c=$c ca=$ca"; end
                if !isapprox(d,da;atol=atold,rtol=rtold); s="$s bad d=$d da=$da"; end
                if s==""; s="ok"; end
                s != "ok" && @info "($a,$b,$c,$d) => ($an,$bn,$cn,$dn) => ($aa,$ba,$ca,$da) : $s"
            end
        end
    # pfrom,pto,α,β,γ,δ)
    end
end
# checkstabledists_paramconv()

function checkmode()
    maxdiff=0.0; maxdiffr=0.0; id=(0,0); idr=(0,0)
    for a in 0.05:0.05:2.0
        for b in 0.0:0.05:1.0
            m=mode(Stable(a,b))
            mopt=mmo1(a,b)
            diff=mopt-m
            diffr=(mopt-m)/max(m,1.e-9)
            if  diff>maxdiff;  maxdiff=diff;   id=(a,b); end
            if diffr>maxdiffr; maxdiffr=diffr; idr=(a,b); end
            if diff > 0.0001
                @info "a=$a b=$b m=$m mopt=$mopt diff=$diff diffr=$diffr"
            end
            if a == 2.0; break; end
        end
    end
    @info "id=$id maxdiff=$maxdiff idr=$idr maxdiffr=$maxdiffr"
end
# checkmode()


"""
    plotnolan_fig1_2 : plot Nolan' fig 1.2, "Stable Densities S(α,0.5;0), parametrization 0, α=0.5,0.75,1,1.25,1.5"
      see figure 1.2 of Chapter 1 (p.8) of Univariate Stable Distributions - John P. Nolan (2020)
"""
function plot_nolan_fig1_2()
    p=plot(title="Stable Densities S(α,0.5) parametrization 0",xlims=(-5.0,5.0),xlabel="α")
    β=0.5
    for α in (0.5,0.75,1.0,1.25,1.5)
        (a,b,c,d)=stabledists_paramconv(0,1,α,β,1.0,0.0)
        d=Stable(a,b,c,d)
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
        d=Stable(α,β)
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
        (a,b,c,d)=stabledists_paramconv(2,1,α,β,1.0,0.0)
        d=Stable(a,b,c,d)
        p=plot!(p,x->pdf(d,x),label="α=$α")
    end
    display(p)
    return nothing
end

function pb0()
    (a,b,c,d)=stabledists_paramconv(0,1,0.75,0.5,1.0,0.0)
    # (0.75, 0.5, 1.0, -1.2071067811865475)
    d=Stable(a,b,c,d)
    plot(x->pdf(d,x))
end
function pb2()
    (a,b,c,d)=stabledists_paramconv(2,1,0.75,0.5,1.0,0.0)
    # (0.75, 0.5, 1.467523221730945, -1.2459371667983017)
    d=Stable(a,b,c,d)
    plot(x->pdf(d,x))
end

function check_mmv2vsv1()
    ar=0.05:0.05:2.0
    br=0.0:0.025:1.0
    for cubic in (false,true)
        maxdiff=0.0; maxloc=(0,0)
        for a in ar
            for b in br
                vmmv1=mmv1(a,b)
                vmmv2=mmv2(a,b;cubic=cubic)
                diff=vmmv2-vmmv1
                if diff>maxdiff
                    maxdiff=diff
                    maxloc=(a,b)
                end
            end
        end
        @info "check_mmv2vsv1 : cubic=$cubic maxdiff=$maxdiff maxloc=$maxloc"
    end
end
# check_mmv2vsv1()

function check_mmv3vsv1()
    ar=0.05:0.05:2.0
    br=0.0:0.025:1.0
    maxdiff=0.0; maxloc=(0,0)
    for a in ar
        for b in br
            vmmv1=mmv1(a,b)
            vmmv3=mmv3(a,b)
            diff=vmmv3-vmmv1
            if diff>maxdiff
                maxdiff=diff
                maxloc=(a,b)
            end
        end
    end
    @info "check_mmv3vsv1 : maxdiff=$maxdiff maxloc=$maxloc"
end
# check_mmv3vsv1()


nothing
