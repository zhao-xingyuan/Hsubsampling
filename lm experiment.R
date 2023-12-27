# RUB: MUST.OW
RUB = function(x,b, m)
{
  xb = sample(x, b, replace=F)
  xrub = sample(xb,m, replace=T)
}
WO = function(x,b, m)
{
  xb = sample(x, b, replace=T)
  xWO = sample(xb,m, replace=F)
}
WW = function(x,b, m)
{
  xb = sample(x, b, replace=T)
  xWW = sample(xb,m, replace=T)
}
WR  =  function(x, m)  sample(x, m, replace=T)
WOR =  function(x, m)  sample(x, m, replace=F)
Poisson = function(x,q)
{
  xpoisson=c()
  for (s in 1:length(x)){
    A=runif(1)
    if (A<q) xpoisson=c(xpoisson,x[s])
  }
  return(xpoisson)
}

eps_prime=0.01
b=200
m=100
n=1000
C=3
eps_WR=log(1+(exp(eps_prime)-1)/(1-(1-1/n)^m))
eps_WOR=log(1+(exp(eps_prime)-1)/(m/n))
eps_Poisson=eps_WOR
sigma <- function(sensitivity, eps){
  sqrt(2*log(1.25/(1/n)))*sensitivity/(eps*n)
}
sigma_WR=sigma(C, eps_WR)
sigma_WOR=sigma(C, eps_WOR)
sigma_Poisson=sigma_WOR

theta_WOR=C/sigma_WOR
theta_WR=C/sigma_WR
theta_Poisson=C/sigma_Poisson

delta_Gau_WB <- function(eps,theta) {
  del=pnorm(theta/2-eps/theta)-exp(eps)*pnorm(-theta/2-eps/theta)
  return (min(del,1))
}
# H_mean_delta=H_var_delta, WW_mean_delta=WW_var_delta
# sensitivity cancelled out in the theta formula, per-query delta' is irrelevant with statistic/sensitivity
# only one per-query delta' for each subsampling method
eps_H=log(1+(exp(eps_prime)-1)/((b/n)*(1-(1-1/b)^m)))
sigma_H=sigma(C, eps_H)
theta_H=C/sigma_H
terms=rep(0,m)
for (k in 1:m){
    delta_k=delta_Gau_WB(eps_H,k*theta_H)
    terms[k]=(b/n)*choose(m,k)*((1/b)^k)*((1-1/b)^(m-k))*delta_k
  }
H_delta=sum(terms)


j=seq(1,b,1)
eta_WW=sum(choose(b,j)*(1/n)^j*(1-1/n)^(b-j)*(1-(1-j/b)^m))
eps_WW=log(1+(exp(eps_prime)-1)/eta_WW)
sigma_WW=sigma(C, eps_WW)
theta_WW=C/sigma_WW
terms=rep(0,b)
for (k in 1:b){
    uterms=rep(0,m)
    for (u in 1:m){
      delta_u=delta_Gau_WB(eps_WW,u*theta_WW)
      uterms[u]=choose(m,u)*(k/b)^u*(1-k/b)^(m-u)*delta_u
    }
    terms[k]=choose(b,k)*(1/n)^k*(1-1/n)^(b-k)*sum(uterms)
  }
WW_delta=sum(terms)

terms=rep(0,m)
for (k in 1:m){
  delta_k=delta_Gau_WB(eps_WR,k*theta_WR)
  terms[k]=choose(m,k)*(1/n)^k*(1-1/n)^(m-k)*delta_k
}
WR_delta=sum(terms)


WOR_delta=(m/n)*(1/n)
Poisson_delta=(m/n)*(1/n)


sgd <- function(
  par,  stepsize = 0.05,
  subsample ='rub', C=1, TT= 100, eps_prime=0.1, n = 1000, b= 200, m= 100)
{
  n  = 1000
  x1 = rnorm(n)
  x2 = rnorm(n)
  y  = 1 + .5*x1 + .2*x2 + rnorm(n)
  X  = cbind(Intercept = 1, x1, x2)
  #lm(y~x1+x2)
  #  1.0300          0.5177       0.1631 
  
  bb = par
  loss = rep(0,TT)                                    
  compute_sigma <- function(sensitivity, eps){
    sqrt(2*log(1.25/(1/n)))*sensitivity/(eps*n)
  }
  
  for (t in 1:TT) {
    if(subsample=='rub') {
      s= RUB(1:n,b,m)
      eps_H=log(1+(exp(eps_prime)-1)/((b/n)*(1-(1-1/b)^m)))
      sigma=compute_sigma(C, eps_H)
    }
    if(subsample=='WO') {
      s= WO(1:n,b,m)
      WO_terms=rep(0,b)
      for (j in 1:b){
        WO_terms[j]=choose(b,j)*(1/n)^j*(1-1/n)^(b-j)*(1-choose(b-j,m)/choose(b,m))
      }
      eta_WO=sum(WO_terms)
      eps_WO=log(1+(exp(eps_prime)-1)/eta_WO)
      sigma=compute_sigma(C, eps_WO)
    }
    if(subsample=='WW') {
      s= RUB(1:n,b,m)
      WW_terms=rep(0,b)
      for (j in 1:b){
        WW_terms[j]=choose(b,j)*(1/n)^j*(1-1/n)^(b-j)*(1-(1-j/b)^m)
      }
      eta_WW=sum(WW_terms)
      eps_WW=log(1+(exp(eps_prime)-1)/eta_WW)
      sigma=compute_sigma(C, eps_WW)
    }
    if(subsample=='WR') {
      s= WR(1:n,m)
      eps_WR=log(1+(exp(eps_prime)-1)/(1-(1-1/n)^m))
      sigma=compute_sigma(C, eps_WR)
    }
    if(subsample=='WOR') {
      s= WOR(1:n,m)  
      eps_WOR=log(1+(exp(eps_prime)-1)/(m/n))
      sigma=compute_sigma(C, eps_WOR)
    }
    # non-private baseline
    if(subsample=='Poisson') {
      s= Poisson(1:n,m/n) 
      sigma=0
    }
    if(subsample=='DP-Poisson') {
      s= Poisson(1:n,m/n) 
      eps_Poisson=log(1+(exp(eps_prime)-1)/(m/n))
      sigma=compute_sigma(C, eps_Poisson)
    }
    grad= matrix(0,length(s),3)
    
    for(i in 1:length(s)){
      Xi   = X[s[i], ]; 
      yi   = y[s[i]]; 
      LP   = Xi %*%bb  
      gradi = 2*Xi %*% (LP - yi); 
      if (subsample=='Poisson') {
        grad[i,] = t(gradi)
      }
      else{grad[i,] = t(gradi/max(1,norm(gradi, type="2")/C))}
    }
    grads=rep(0,3)
    for(k in 1:3){
        grads[k] = sum(grad[,k])/length(s)+rnorm(1,0, sigma)
    } 
    bb = bb - stepsize*grads
    loss[t]     = sum((X %*%bb - y)^2)/n
  }
  list(est= bb,loss=loss, sig=sigma, rmse=sqrt(mean(loss[(TT-15):TT])))
}


set.seed(1234)
n  = 1000
x1 = rnorm(n)
x2 = rnorm(n)
y  = 1 + .5*x1 + .2*x2 + rnorm(n)
X  = cbind(Intercept = 1, x1, x2)
lm(y~x1+x2)
#  1.0300          0.5177       0.1631  

# training loss in a single dataset
TT=200
set.seed(9)
x11()
par(mfrow=c(1,2), mar=c(4,4,1,1))

out=sgd(rep(-1,3), subsample="WR", m=100, stepsize=0.04, TT=TT, C=3,eps_prime=0.01)
out$sig;plot(1:TT,out$loss,type="l", xlab="", ylab="", ylim=c(0,9),col='blue',lty=2, lwd = 2) 
out=sgd(rep(-1,3), subsample="WOR", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.01)
out$sig; lines(1:TT,out$loss,type="l", col='green3',lty='aa', lwd = 2)
out=sgd(rep(-1,3), subsample="rub", m=100, stepsize=0.04, TT=TT, C=3,eps_prime=0.01)
out$sig; lines(1:TT,out$loss,type="l", col='red',lty=3, lwd = 2)
out=sgd(rep(-1,3), subsample="WW", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.01)
out$sig; lines(1:TT,out$loss,type="l", col='purple',lty=5, lwd = 2)
out=sgd(rep(-1,3), subsample="DP-Poisson", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.01)
out$sig; lines(1:TT,out$loss,type="l", col='brown',lty=6, lwd = 2)
out=sgd(rep(-1,3), subsample="Poisson", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.01)
lines(1:TT,out$loss,type="l", col='black',lty=1, lwd = 2)
legend(15, 9, c("Poisson (non-private)","Poisson","WOR", "WR", "MUST.OW","MUST.WW"),  col=c('black','brown','green3','blue','red','purple'),lty=c('solid','twodash','aa','dashed','dotted','longdash'),bty = "n",cex=1.1)
mtext("loss",2,2)
mtext("iteration",1,2,cex=1)

out=sgd(rep(-1,3), subsample="WR", m=100, stepsize=0.04, TT=TT, C=3,eps_prime=0.001)
out$sig; plot(1:TT,out$loss,type="l", xlab="", ylab="", ylim=c(0,9),col='blue',lty=2, lwd = 2) 
out=sgd(rep(-1,3), subsample="WOR", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.001)
out$sig; lines(1:TT,out$loss,type="l", col='green3',lty='aa', lwd = 2)
out=sgd(rep(-1,3), subsample="rub", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.001)
out$sig; lines(1:TT,out$loss,type="l", col='red',lty=3, lwd = 2)
out=sgd(rep(-1,3), subsample="WW", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.001)
out$sig; lines(1:TT,out$loss,type="l", col='purple',lty=5, lwd = 2)
out=sgd(rep(-1,3), subsample="DP-Poisson", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.001)
out$sig; lines(1:TT,out$loss,type="l", col='brown',lty=6, lwd = 2)
out=sgd(rep(-1,3), subsample="Poisson", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.001)
lines(1:TT,out$loss,type="l", col='black',lty=1, lwd = 2)
mtext("loss",2,2)
mtext("iteration",1,2,cex=1)


rmse=matrix(0,13,2)
num_repeats=200
rmse1=rep(0,num_repeats)
rmse05=matrix(0,num_repeats,6)
rmse01=matrix(0,num_repeats,6)

true_bb=c(1,0.5,0.2)
beta_bias=array(0,c(num_repeats,13,3))
set.seed(9)
start_time<-Sys.time()
for(i in 1:num_repeats){
  ### eps_prime=0.1
  out=sgd(rep(-1,3), subsample="Poisson", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.1)
  rmse1[i]= out$rmse;beta_bias[i,1,]=out$est-true_bb
  out=sgd(rep(-1,3), subsample="WR", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.01)
  rmse05[i,1]= out$rmse;beta_bias[i,4,]=out$est-true_bb
  out=sgd(rep(-1,3), subsample="WOR", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.01) 
  rmse05[i,2]= out$rmse; beta_bias[i,3,]=out$est-true_bb
  out=sgd(rep(-1,3), subsample="DP-Poisson", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.01) 
  rmse05[i,3]= out$rmse; beta_bias[i,2,]=out$est-true_bb
  out=sgd(rep(-1,3), subsample="rub", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.01)
  rmse05[i,4]= out$rmse;beta_bias[i,6,]=out$est-true_bb
  out=sgd(rep(-1,3), subsample="WO", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.01) 
  rmse05[i,5]= out$rmse; beta_bias[i,5,]=out$est-true_bb
  out=sgd(rep(-1,3), subsample="WW", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.01)
  rmse05[i,6]= out$rmse;beta_bias[i,7,]=out$est-true_bb
  ### eps_prime=0.001
  out=sgd(rep(-1,3), subsample="WR", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.001)
  rmse01[i,1]= out$rmse;beta_bias[i,10,]=out$est-true_bb
  out=sgd(rep(-1,3), subsample="WOR", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.001) 
  rmse01[i,2]= out$rmse;beta_bias[i,9,]=out$est-true_bb
  out=sgd(rep(-1,3), subsample="DP-Poisson", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.001) 
  rmse01[i,3]= out$rmse; beta_bias[i,8,]=out$est-true_bb
  out=sgd(rep(-1,3), subsample="rub", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.001)
  rmse01[i,4]= out$rmse;beta_bias[i,12,]=out$est-true_bb
  out=sgd(rep(-1,3), subsample="WO", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.001) 
  rmse01[i,5]= out$rmse; beta_bias[i,11,]=out$est-true_bb
  out=sgd(rep(-1,3), subsample="WW", m=100, stepsize=0.04, TT=TT, C=3, eps_prime=0.001)
  rmse01[i,6]= out$rmse;beta_bias[i,13,]=out$est-true_bb

}
end_time<-Sys.time()
end_time-start_time
#Time difference of 1.183985 hours

rmse[1,1]=mean(rmse1) 
rmse[1,2]=sd(rmse1)
rmse[2:7,1]=apply(rmse05,2,mean) 
rmse[2:7,2]=apply(rmse05, 2, sd)
rmse[8:13,1]=apply(rmse01,2,mean) 
rmse[8:13,2]=apply(rmse01, 2, sd)

beta_bias_lm=apply(beta_bias,c(2,3),mean)
beta_rmse_lm=sqrt(apply(beta_bias^2,c(2,3),mean))


write.csv(rmse, "utility lm rmse_T200.csv")
write.csv(beta_bias_lm, "utility lm beta bias_T200.csv")
write.csv(beta_rmse_lm, "utility lm beta rmse_T200.csv")
