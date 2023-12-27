eps_prime=0.1
repeats=200
B= 500
n= 300
m= 30
bound = 4
b0= c(10, 20, 30, 50, 70, 100); L=length(b0)
WOR=rep(0,repeats)
WR=rep(0,repeats)
WOR_var=rep(0,repeats)
WR_var=rep(0,repeats)
Poisson_mean=rep(0,repeats)
Poisson_var=rep(0,repeats)
nonPrivate_mean=rep(0,repeats)
nonPrivate_var=rep(0,repeats)
RUB_mean=matrix(0,L,repeats)
RUB_var=matrix(0,L,repeats)
WW_mean=matrix(0,L,repeats)
WW_var=matrix(0,L,repeats)



eps_WR=log(1+(exp(eps_prime)-1)/(1-(1-1/n)^m))
eps_WOR=log(1+(exp(eps_prime)-1)/(m/n))
eps_Poisson=eps_WOR
sigma <- function(sensitivity, eps){
  sqrt(2*log(1.25/(1/n)))*sensitivity/eps
}
sigma_WR=sigma(2*bound/n, eps_WR)
sigma_WOR=sigma(2*bound/n, eps_WOR)
sigma_WR_var=sigma((2*bound)^2/n, eps_WR)
sigma_WOR_var=sigma((2*bound)^2/n, eps_WOR)
sigma_Poisson=sigma_WOR
sigma_Poisson_var=sigma_WOR_var

theta_WOR=(2*bound/n)/sigma_WOR
theta_WOR_var=((2*bound)^2/n)/sigma_WOR_var
theta_WR=(2*bound/n)/sigma_WR
theta_WR_var=((2*bound)^2/n)/sigma_WR_var
theta_Poisson=(2*bound/n)/sigma_Poisson
theta_Poisson_var=((2*bound)^2/n)/sigma_Poisson_var

sigma_H=rep(0,L)
sigma_WW=rep(0,L)
sigma_H_var=rep(0,L)
sigma_WW_var =rep(0,L)  
eps_H=rep(0,L)
eps_WW =rep(0,L)
theta_H=rep(0,L)
theta_H_var =rep(0,L)
H_mean_delta=rep(0,L)
H_var_delta=rep(0,L)
theta_WW=rep(0,L)
theta_WW_var =rep(0,L)
WW_mean_delta=rep(0,L)
WW_var_delta=rep(0,L)

delta_Gau_WB <- function(eps,theta) {
  del=pnorm(theta/2-eps/theta)-exp(eps)*pnorm(-theta/2-eps/theta)
  return (min(del,1))
}
# H_mean_delta=H_var_delta, WW_mean_delta=WW_var_delta
# sensitivity cancelled out in the theta formula, per-query delta' is irrelevant with statistic/sensitivity
# only one per-query delta' for each subsampling method
for (l in 1:L){
  b= b0[l]
  eps_H[l]=log(1+(exp(eps_prime)-1)/((b/n)*(1-(1-1/b)^m)))
  sigma_H[l]=sigma(2*bound/n, eps_H[l])
  sigma_H_var[l]=sigma((2*bound)^2/n, eps_H[l])
  theta_H[l]=(2*bound/n)/sigma_H[l]
  theta_H_var[l]=((2*bound)^2/n)/sigma_H_var[l]
  terms=rep(0,m)
  for (k in 1:m){
    delta_k=delta_Gau_WB(eps_H[l],k*theta_H[l])
    terms[k]=(b/n)*choose(m,k)*((1/b)^k)*((1-1/b)^(m-k))*delta_k
  }
  H_mean_delta[l]=sum(terms)
  var_terms=rep(0,m)
  for (k in 1:m){
    delta_k=delta_Gau_WB(eps_H[l],k*theta_H_var[l])
    var_terms[k]=(b/n)*choose(m,k)*((1/b)^k)*((1-1/b)^(m-k))*delta_k
  }
  H_var_delta[l]=sum(var_terms)
}

for (l in 1:L){
  b= b0[l]
  j=seq(1,b,1)
  eta_WW=sum(choose(b,j)*(1/n)^j*(1-1/n)^(b-j)*(1-(1-j/b)^m))
  eps_WW[l]=log(1+(exp(eps_prime)-1)/eta_WW)
  sigma_WW[l]=sigma(2*bound/n, eps_WW[l])
  sigma_WW_var[l]=sigma((2*bound)^2/n, eps_WW[l])
  theta_WW[l]=(2*bound/n)/sigma_WW[l]
  theta_WW_var[l]=((2*bound)^2/n)/sigma_WW_var[l]
  terms=rep(0,b)
  for (k in 1:b){
    uterms=rep(0,m)
    for (u in 1:m){
      delta_u=delta_Gau_WB(eps_WW[l],u*theta_WW[l])
      uterms[u]=choose(m,u)*(k/b)^u*(1-k/b)^(m-u)*delta_u
    }
    terms[k]=choose(b,k)*(1/n)^k*(1-1/n)^(b-k)*sum(uterms)
  }
  WW_mean_delta[l]=sum(terms)
  var_terms=rep(0,b)
  for (k in 1:b){
    uterms=rep(0,m)
    for (u in 1:m){
      delta_u=delta_Gau_WB(eps_WW[l],u*theta_WW_var[l])
      uterms[u]=choose(m,u)*(k/b)^u*(1-k/b)^(m-u)*delta_u
    }
    var_terms[k]=choose(b,k)*(1/n)^k*(1-1/n)^(b-k)*sum(uterms)
  }
  WW_var_delta[l]=sum(var_terms)
}

terms=rep(0,m)
for (k in 1:m){
  delta_k=delta_Gau_WB(eps_WR,k*theta_WR)
  terms[k]=choose(m,k)*(1/n)^k*(1-1/n)^(m-k)*delta_k
  }
WR_mean_delta=sum(terms)

terms=rep(0,m)
for (k in 1:m){
  delta_k=delta_Gau_WB(eps_WR,k*theta_WR_var)
  terms[k]=choose(m,k)*(1/n)^k*(1-1/n)^(m-k)*delta_k
}
WR_var_delta=sum(terms)


WOR_mean_delta=(m/n)*(1/n)

Poisson_mean_delta=(m/n)*(1/n)

Poisson = function(x,q)
{
  xpoisson=c()
  for (s in 1:length(x)){
    A=runif(1)
    if (A<q) xpoisson=c(xpoisson,x[s])
  }
  return(xpoisson)
}

set.seed(10)
start_time <- Sys.time()
for (t in 1:repeats){
  x= rnorm(n, 0,1)
  
  estWOR =rep(0, B)
  varWOR =rep(0, B)
  for(i in 1:B)
    {
      xb = sample(x,m, replace=F)
      estWOR[i] = mean(xb)+rnorm(1,0,sigma_WOR)
      varWOR[i] = var(xb)+rnorm(1,0,sigma_WOR_var)
    }
    WOR[t]=mean(estWOR)
    WOR_var[t]=mean(varWOR)
    
    estPoisson =rep(0, B)
    varPoisson =rep(0, B)
    est_nonprivate =rep(0, B)
    var_nonprivate =rep(0, B)
    for(i in 1:B)
    {
      xb = Poisson(x,m/n)
      est_nonprivate[i] = mean(xb)
      var_nonprivate[i] = var(xb)
      estPoisson[i] = mean(xb)+rnorm(1,0,sigma_Poisson)
      varPoisson[i] = var(xb)+rnorm(1,0,sigma_Poisson_var)
    }
    Poisson_mean[t]=mean(estPoisson)
    Poisson_var[t]=mean(varPoisson)
    nonPrivate_mean[t]=mean(est_nonprivate)
    nonPrivate_var[t]=mean(var_nonprivate)
    
    estWR =rep(0, B) # subsample bootstrap
    varWR =rep(0, B)
    for(i in 1:B)
    {
      xb = sample(x, m, replace=T)
      estWR[i] = mean(xb)+rnorm(1,0,sigma_WR)
      varWR[i] = var(xb)+rnorm(1,0,sigma_WR_var)
    }
    WR[t]=mean(estWR)
    WR_var[t]=mean(varWR)
  
    sigma_H=rep(0,L)
    sigma_WW=rep(0,L)
    sigma_H_var=rep(0,L)
    sigma_WW_var =rep(0,L)  
    eps_H=rep(0,L)
    eps_WW =rep(0,L)
  for (l in 1:L){
    b= b0[l]
    eps_H[l]=log(1+(exp(eps_prime)-1)/((b/n)*(1-(1-1/b)^m)))
    sigma_H[l]=sigma(2*bound/n, eps_H[l])
    sigma_H_var[l]=sigma((2*bound)^2/n, eps_H[l])
    estRUB =rep(0, B)
    varRUB =rep(0, B)
    for(i in 1:B)
    {
      xb = sample(x, b, replace=F)
      xrub = sample(xb,m, replace=T)
      estRUB[i] = mean(xrub)+rnorm(1,0,sigma_H[l]) 
      varRUB[i] = var(xrub)+rnorm(1,0,sigma_H_var[l]) 
    }
    RUB_mean[l,t]=mean(estRUB)
    RUB_var[l,t]=mean(varRUB)
  
    WW_terms=rep(0,b)
    for (j in 1:b){
        WW_terms[j]=choose(b,j)*(1/n)^j*(1-1/n)^(b-j)*(1-(1-j/b)^m)
      }
    eta_WW=sum(WW_terms)
    eps_WW[l]=log(1+(exp(eps_prime)-1)/eta_WW)
    sigma_WW[l]=sigma(2*bound/n, eps_WW[l])
    sigma_WW_var[l]=sigma((2*bound)^2/n, eps_WW[l])
    estWW =rep(0, B)
    varWW =rep(0, B)
    for(i in 1:B)
      {
        xb = sample(x, b, replace=T)
        xWW = sample(xb,m, replace=T)
        estWW[i] = mean(xWW)+rnorm(1,0,sigma_WW[l]) 
        varWW[i] = var(xWW)+rnorm(1,0,sigma_WW_var[l]) 
    }
    WW_mean[l,t]=mean(estWW)
    WW_var[l,t]=mean(varWW)
    }
}
end_time <- Sys.time()
end_time - start_time

# point estimate and SD for sample mean and sample variance
mean(WOR)
sd(WOR)
mean(WOR_var)
sd(WOR_var)

mean(WR)
sd(WR)
mean(WR_var)
sd(WR_var)

mean(Poisson_mean)
sd(Poisson_mean)
mean(Poisson_var)
sd(Poisson_var)

mean(nonPrivate_mean)
sd(nonPrivate_mean)
mean(nonPrivate_var)
sd(nonPrivate_var)

apply(RUB_mean,1,mean)
apply(RUB_var,1,mean)
apply(RUB_mean,1,sd)
apply(RUB_var,1,sd)

apply(WW_mean,1,mean)
apply(WW_var,1,mean)
apply(WW_mean,1,sd)
apply(WW_var,1,sd)




