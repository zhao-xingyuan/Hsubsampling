library(latex2exp)

# (a) epsilon PA vs M plot
n=1000
m=400
b=500


eps=seq(0.01,6,0.01)
eps_WOR=log(1+(m/n)*(exp(eps)-1))
eps_WR=log(1+(1-(1-1/n)^m)*(exp(eps)-1))
eps_H=log(1+((b/n)*(1-(1-1/b)^m))*(exp(eps)-1))
j=seq(1,b,1)
eta_WO=sum(choose(b,j)*(1/n)^j*(1-1/n)^(b-j)*(1-choose(b-j,m)/choose(b,m)))
eps_WO=log(1+eta_WO*(exp(eps)-1))
eta_WW=sum(choose(b,j)*(1/n)^j*(1-1/n)^(b-j)*(1-(1-j/b)^m))
eps_WW=log(1+eta_WW*(exp(eps)-1))
x11()
plot(eps,eps,pch='',xlim=c(0,6),ylim=c(0,6),xlab=expression(epsilon),ylab=expression(paste(epsilon,"'")),cex.lab = 1.5)
lines(eps,eps, col="black",lty=1,lwd=2) # no subsampling
lines(eps,eps_WOR, col="green3",lty=2,lwd=2)
lines(eps,eps_WR, col="blue",lty=3,lwd=2)
lines(eps,eps_H, col="red",lty=4,lwd=2)
lines(eps,eps_WW, col="purple",lty=6,lwd=2)
legend(0,6,legend=c(TeX("$\\textit{M}$"),TeX("$\\textit{M}$+WOR"),TeX("$\\textit{M}$+WR"),TeX("$\\textit{M}$+MUST.OW"),TeX("$\\textit{M}$+MUST.WW")), col=c("black","green3","blue","red", "purple"),
       bty = 'n' ,lwd=c(2,2,2,2,2),lty=c(1,2,3,4,6),ncol=1,cex=1.2)

# (b) delta PA vs M plot - Lap
theta=0.25
delta = ifelse(1-exp((eps-theta)/2)<0,0,1-exp((eps-theta)/2))
delta_Lap <- function(eps) {
  ifelse(1-exp((eps-theta)/2)<0,0,1-exp((eps-theta)/2))
}
delta_Lap_WB <- function(eps,theta) {
  ifelse(1-exp((eps-theta)/2)<0,0,1-exp((eps-theta)/2))
}
eps_WOR=log(1+(m/n)*(exp(eps)-1))
eps_WR=log(1+(1-(1-1/n)^m)*(exp(eps)-1))
eps_H=log(1+((b/n)*(1-(1-1/b)^m))*(exp(eps)-1))
j=seq(1,b,1)
eta_WO=sum(choose(b,j)*(1/n)^j*(1-1/n)^(b-j)*(1-choose(b-j,m)/choose(b,m)))
eps_WO=log(1+eta_WO*(exp(eps)-1))
eta_WW=sum(choose(b,j)*(1/n)^j*(1-1/n)^(b-j)*(1-(1-j/b)^m))
eps_WW=log(1+eta_WW*(exp(eps)-1))

delta_WR_WB=rep(0,length(eps_WR))
for (i in 1:length(eps_WR)){
  terms=rep(0,m)
  for (k in 1:m){
    delta_k=delta_Lap_WB(eps[i],k*theta)
    terms[k]=choose(m,k)*(1/n)^k*(1-1/n)^(m-k)*delta_k
  }
  delta_WR_WB[i]=sum(terms)
}

delta_H_WB=rep(0,length(eps_H))
for (i in 1:length(eps_H)){
  terms=rep(0,m)
  for (k in 1:m){
    delta_k=delta_Lap_WB(eps[i],k*theta)
    terms[k]=(b/n)*choose(m,k)*((1/b)^k)*((1-1/b)^(m-k))*delta_k
  }
  delta_H_WB[i]=sum(terms)
}
delta_WOR_WB=rep(0,length(eps_WOR))
for (i in 1:length(eps_WOR)){
  delta_WOR_WB[i]=(m/n)*delta_Lap_WB(eps[i],theta)
}

delta_WO_WB=rep(0,length(eps_WO))
for (i in 1:length(eps_WO)){
  terms=rep(0,b)
  for (k in 1:m){
    uterms=rep(0,k)
    for (u in 1:k){
      delta_u=delta_Lap_WB(eps[i],u*theta)
      uterms[u]=(choose(k,u)*choose(b-k,m-u)/choose(b,m))*delta_u
    }
    terms[k]=choose(b,k)*(1/n)^k*(1-1/n)^(b-k)*sum(uterms)
  }
  for (k in (m+1):b){
    umterms=rep(0,m)
    for (u in 1:m){
      delta_u=delta_Lap_WB(eps[i],u*theta)
      umterms[u]=(choose(k,u)*choose(b-k,m-u)/choose(b,m))*delta_u
    }
    terms[k]=choose(b,k)*(1/n)^k*(1-1/n)^(b-k)*sum(umterms)
  }
  delta_WO_WB[i]=sum(terms)
}

delta_WW_WB=rep(0,length(eps_WW))
for (i in 1:length(eps_WW)){
  terms=rep(0,b)
  for (k in 1:b){
    uterms=rep(0,m)
    for (u in 1:m){
      delta_u=delta_Lap_WB(eps[i],u*theta)
      uterms[u]=choose(m,u)*(k/b)^u*(1-k/b)^(m-u)*delta_u
    }
    terms[k]=choose(b,k)*(1/n)^k*(1-1/n)^(b-k)*sum(uterms)
  }
  delta_WW_WB[i]=sum(terms)
}

x11()
plot(delta,delta,pch='',xlim=c(0,0.4),ylim=c(0,0.4),xlab=expression(delta),ylab=expression(paste(delta,"'")),cex.lab = 1.5)
lines(delta,delta, col="black",lty=1,lwd=2) # no subsampling
lines(delta,delta_WOR_WB, col="green3",lty=2,lwd=2)
lines(delta,delta_WR_WB, col="blue",lty=3,lwd=2)
lines(delta,delta_H_WB, col="red",lty=4,lwd=2)
lines(delta,delta_WO_WB, col="orange",lty=5,lwd=2)
lines(delta,delta_WW_WB, col="purple",lty=6,lwd=2)
legend(-0.02,0.42,legend=c("Laplace","Laplace+WOR","Laplace+WR","Laplace+MUST.OW","Laplace+MUST.WO","Laplace+MUST.WW"), col=c("black","green3","blue","red","orange", "purple"),
       bty = 'n' ,lwd=c(2,2,2,2,2,2),lty=c(1,2,3,4,5,6),ncol=1,cex=1.2) 

PA_Lap_results=cbind(eps,eps_WOR, eps_WR, eps_H, eps_WO, eps_WW, delta,delta_WOR_WB,delta_WR_WB,delta_H_WB,delta_WO_WB,delta_WW_WB)
write.csv(PA_Lap_results, "PA_Laplace0.25.csv")
# delta/sigma=1
x11()
plot(delta,delta,pch='',xlim=c(0,0.4),ylim=c(0,0.4),xlab=expression(delta),ylab=expression(paste(delta,"'")),cex.lab = 1.5)
lines(delta,delta, col="black",lty=1,lwd=2) # no subsampling
lines(delta,delta_WOR_WB, col="green3",lty=2,lwd=2)
lines(delta,delta_WR_WB, col="blue",lty=3,lwd=2)
lines(delta,delta_H_WB, col="red",lty=4,lwd=2)
lines(delta,delta_WW_WB, col="purple",lty=6,lwd=2)
legend(-0.02,0.42,legend=c("Laplace","Laplace+WOR","Laplace+WR","Laplace+MUST.OW","Laplace+MUST.WW"), col=c("black","green3","blue","red", "purple"),
       bty = 'n' ,lwd=c(2,2,2,2,2),lty=c(1,2,3,4,6),ncol=1,cex=1.2) 
# delta/sigma=0.25
x11()
plot(delta,delta,pch='',xlim=c(0,0.11),ylim=c(0,0.11),xlab=expression(delta),ylab=expression(paste(delta,"'")),cex.lab = 1.5)
lines(delta,delta, col="black",lty=1,lwd=2) # no subsampling
lines(delta,delta_WOR_WB, col="green3",lty=2,lwd=2)
lines(delta,delta_WR_WB, col="blue",lty=3,lwd=2)
lines(delta,delta_H_WB, col="red",lty=4,lwd=2)
lines(delta,delta_WW_WB, col="purple",lty=6,lwd=2)
legend(0,0.11,legend=c("Laplace","Laplace+WOR","Laplace+WR","Laplace+MUST.OW","Laplace+MUST.WW"), col=c("black","green3","blue","red", "purple"),
       bty = 'n' ,lwd=c(2,2,2,2,2),lty=c(1,2,3,4,6),ncol=1,cex=1.2) 

# eps'/eps vs delta'-delta plot 
#theta=1
x11()
plot(eps_WOR/eps,delta_WOR_WB-delta,lty=2,lwd=2,pch='',col="green3",xlim=c(0.2,0.9),ylim=c(-0.3,0.3),xlab=expression(paste(epsilon,"'","/",epsilon)),ylab=expression(paste(delta,"'","-",delta)),cex.lab = 1.5)
lines(eps_WOR/eps,delta_WOR_WB-delta, col="green3",lty=2,lwd=2)
lines(eps_WR/eps,delta_WR_WB-delta, col="blue",lty=3,lwd=2)
lines(eps_H/eps,delta_H_WB-delta, col="red",lty=4,lwd=2)
lines(eps_WW/eps,delta_WW_WB-delta, col="purple",lty=6,lwd=2)
legend(0.2,0.26,legend=c("Laplace+WOR","Laplace+WR","Laplace+MUST.OW","Laplace+MUST.WW"), col=c("green3","blue","red", "purple"),
       bty = 'n' ,lwd=c(2,2,2,2),lty=c(2,3,4,6),ncol=1,cex=1.2) 
# theta=0.25
plot(eps_WOR/eps,delta_WOR_WB-delta,lty=2,lwd=2,pch='',col="green3",xlim=c(0.2,0.9),ylim=c(-0.1,0.1),xlab=expression(paste(epsilon,"'","/",epsilon)),ylab=expression(paste(delta,"'","-",delta)),cex.lab = 1.5)
lines(eps_WOR/eps,delta_WOR_WB-delta, col="green3",lty=2,lwd=2)
lines(eps_WR/eps,delta_WR_WB-delta, col="blue",lty=3,lwd=2)
lines(eps_H/eps,delta_H_WB-delta, col="red",lty=4,lwd=2)
lines(eps_WW/eps,delta_WW_WB-delta, col="purple",lty=6,lwd=2)
legend(0.2,0.09,legend=c("Laplace+WOR","Laplace+WR","Laplace+MUST.OW","Laplace+MUST.WW"), col=c("green3","blue","red", "purple"),
       bty = 'n' ,lwd=c(2,2,2,2),lty=c(2,3,4,6),ncol=1,cex=1.2) 

# (c) delta PA vs M plot - Gaussian
theta=1
delta = pnorm(theta/2-eps/theta)-exp(eps)*pnorm(-theta/2-eps/theta)
delta_Gau <- function(eps) {
  del=pnorm(theta/2-eps/theta)-exp(eps)*pnorm(-theta/2-eps/theta)
  return (min(del,1))
}
delta_Gau_WB <- function(eps,theta) {
  del=pnorm(theta/2-eps/theta)-exp(eps)*pnorm(-theta/2-eps/theta)
  return (min(del,1))
}

delta_WR_WB=rep(0,length(eps_WR))
for (i in 1:length(eps_WR)){
  terms=rep(0,m)
  for (k in 1:m){
    delta_k=delta_Gau_WB(eps[i],k*theta)
    terms[k]=choose(m,k)*(1/n)^k*(1-1/n)^(m-k)*delta_k
  }
  delta_WR_WB[i]=sum(terms)
}

delta_H_WB=rep(0,length(eps_H))
for (i in 1:length(eps_H)){
  terms=rep(0,m)
  for (k in 1:m){
    delta_k=delta_Gau_WB(eps[i],k*theta)
    terms[k]=(b/n)*choose(m,k)*((1/b)^k)*((1-1/b)^(m-k))*delta_k
  }
  delta_H_WB[i]=sum(terms)
}

delta_WOR_WB=rep(0,length(eps_WOR))
for (i in 1:length(eps_WOR)){
  delta_WOR_WB[i]=(m/n)*delta_Gau_WB(eps[i],theta)
}

delta_WO_WB=rep(0,length(eps_WO))
for (i in 1:length(eps_WO)){
  terms=rep(0,b)
  for (k in 1:m){
    uterms=rep(0,k)
    for (u in 1:k){
      delta_u=delta_Gau_WB(eps[i],u*theta)
      uterms[u]=(choose(k,u)*choose(b-k,m-u)/choose(b,m))*delta_u
    }
    terms[k]=choose(b,k)*(1/n)^k*(1-1/n)^(b-k)*sum(uterms)
  }
  for (k in (m+1):b){
    umterms=rep(0,m)
    for (u in 1:m){
      delta_u=delta_Gau_WB(eps[i],u*theta)
      umterms[u]=(choose(k,u)*choose(b-k,m-u)/choose(b,m))*delta_u
    }
    terms[k]=choose(b,k)*(1/n)^k*(1-1/n)^(b-k)*sum(umterms)
  }
  delta_WO_WB[i]=sum(terms)
}

delta_WW_WB=rep(0,length(eps_WW))
for (i in 1:length(eps_WW)){
  terms=rep(0,b)
  for (k in 1:b){
    uterms=rep(0,m)
    for (u in 1:m){
      delta_u=delta_Gau_WB(eps[i],u*theta)
      uterms[u]=choose(m,u)*(k/b)^u*(1-k/b)^(m-u)*delta_u
    }
    terms[k]=choose(b,k)*(1/n)^k*(1-1/n)^(b-k)*sum(uterms)
  }
  delta_WW_WB[i]=sum(terms)
}

#Delta/sigma=1
x11()
plot(delta,delta,pch='',xlim=c(0,0.4),ylim=c(0,0.4),xlab=expression(delta),ylab=expression(paste(delta,"'")),cex.lab = 1.5)
lines(delta,delta, col="black",lty=1,lwd=2) # no subsampling
lines(delta,delta_WOR_WB, col="green3",lty=2,lwd=2)
lines(delta,delta_WR_WB, col="blue",lty=3,lwd=2)
lines(delta,delta_H_WB, col="red",lty=4,lwd=2)
lines(delta,delta_WW_WB, col="purple",lty=6,lwd=2)
legend(-0.02,0.42,legend=c("Gaussian","Gaussian+WOR","Gaussian+WR","Gaussian+MUST.OW","Gaussian+MUST.WW"), col=c("black","green3","blue","red", "purple"),
       bty = 'n' ,lwd=c(2,2,2,2,2),lty=c(1,2,3,4,6),ncol=1,cex=1.2) 

PA_Gauss_results=cbind(eps,eps_WOR, eps_WR, eps_H, eps_WO, eps_WW, delta,delta_WOR_WB,delta_WR_WB,delta_H_WB,delta_WO_WB,delta_WW_WB)
write.csv(PA_Gauss_results, "PA_Gaussian0.25.csv")
#Delta/sigma=0.25
x11()
plot(delta,delta,pch='',xlim=c(0,0.09),ylim=c(0,0.09),xlab=expression(delta),ylab=expression(paste(delta,"'")),cex.lab = 1.5)
lines(delta,delta, col="black",lty=1,lwd=2) # no subsampling
lines(delta,delta_WOR_WB, col="green3",lty=2,lwd=2)
lines(delta,delta_WR_WB, col="blue",lty=3,lwd=2)
lines(delta,delta_H_WB, col="red",lty=4,lwd=2)
lines(delta,delta_WW_WB, col="purple",lty=6,lwd=2)
legend(0,0.092,legend=c("Gaussian","Gaussian+WOR","Gaussian+WR","Gaussian+MUST.OW","Gaussian+MUST.WW"), col=c("black","green3","blue","red", "purple"),
       bty = 'n' ,lwd=c(2,2,2,2,2),lty=c(1,2,3,4,6),ncol=1,cex=1.2) 

# eps'/eps vs delta'-delta plot 
#theta=1
x11()
plot(eps_WOR/eps,delta_WOR_WB-delta,lty=2,lwd=2,pch='',col="green3",xlim=c(0.2,0.9),ylim=c(-0.3,0.3),xlab=expression(paste(epsilon,"'","/",epsilon)),ylab=expression(paste(delta,"'","-",delta)),cex.lab = 1.5)
lines(eps_WOR/eps,delta_WOR_WB-delta, col="green3",lty=2,lwd=2)
lines(eps_WR/eps,delta_WR_WB-delta, col="blue",lty=3,lwd=2)
lines(eps_H/eps,delta_H_WB-delta, col="red",lty=4,lwd=2)
lines(eps_WW/eps,delta_WW_WB-delta, col="purple",lty=6,lwd=2)
legend(0.2,0.26,legend=c("Gaussian+WOR","Gaussian+WR","Gaussian+MUST.OW","Gaussian+MUST.WW"), col=c("green3","blue","red", "purple"),
       bty = 'n' ,lwd=c(2,2,2,2),lty=c(2,3,4,6),ncol=1,cex=1.2) 
#theta=0.25
x11()
plot(eps_WOR/eps,delta_WOR_WB-delta,lty=2,lwd=2,pch='',col="green3",xlim=c(0.2,0.9),ylim=c(-0.1,0.1),xlab=expression(paste(epsilon,"'","/",epsilon)),ylab=expression(paste(delta,"'","-",delta)),cex.lab = 1.5)
lines(eps_WOR/eps,delta_WOR_WB-delta, col="green3",lty=2,lwd=2)
lines(eps_WR/eps,delta_WR_WB-delta, col="blue",lty=3,lwd=2)
lines(eps_H/eps,delta_H_WB-delta, col="red",lty=4,lwd=2)
lines(eps_WW/eps,delta_WW_WB-delta, col="purple",lty=6,lwd=2)
legend(0.2,0.08,legend=c("Gaussian+WOR","Gaussian+WR","Gaussian+MUST.OW","Gaussian+MUST.WW"), col=c("green3","blue","red", "purple"),
       bty = 'n' ,lwd=c(2,2,2,2),lty=c(2,3,4,6),ncol=1,cex=1.2) 
