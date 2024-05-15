using RCall
R"library(stabledist)"


function r_dstable(x,alpha,beta,gamma=1.0,delta=0.0,pm=0)
    rcopy(rcall(:dstable,x,alpha,beta,gamma,delta,pm))
end
function r_pstable(q,alpha,beta,gamma=1.0,delta=0.0,pm=0)
    rcopy(rcall(:pstable,q,alpha,beta,gamma,delta,pm))
end
function r_qstable(p,alpha,beta,gamma=1.0,delta=0.0,pm=0)
    rcopy(rcall(:qstable,p,alpha,beta,gamma,delta,pm))
end
function r_rstable(n,alpha,beta,gamma=1.0,delta=0.0,pm=0)
    rcopy(rcall(:rstable,n,alpha,beta,gamma,delta,pm))
end
function r_stableMode(alpha,beta) # gamma = 1, delta = 0, pm = 0
    rcopy(rcall(:stableMode,alpha,beta))
end

function r_dstabler(x,alpha,beta,gamma=1.0,delta=0.0,pm=0)
    redirect_stdio(;stderr=devnull) do
        r_dstable(x,alpha,beta,gamma,delta,pm)
    end
end
function r_pstabler(x,alpha,beta,gamma=1.0,delta=0.0,pm=0)
    redirect_stdio(;stderr=devnull) do
        r_pstable(x,alpha,beta,gamma,delta,pm)
    end
end


nothing
