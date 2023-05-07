#### Figure 1
# (a) epsilon PA vs M plot
n=1000
m=400
b=m^(0.9)

eps=seq(0.01,6,0.01)
eps_WOR=log(1+(m/n)*(exp(eps)-1))
eps_WR=log(1+(1-(1-1/n)^m)*(exp(eps)-1))
eps_H=log(1+((b/n)*(1-(1-1/b)^m))*(exp(eps)-1))
x11()
plot(eps,eps,pch='',xlim=c(0,6),ylim=c(0,6),xlab=expression(epsilon),ylab=expression(paste(epsilon,"'")),cex.lab = 1.5)
lines(eps,eps, col="black",lty=1,lwd=2) # no subsampling
lines(eps,eps_WOR, col="green3",lty=2,lwd=2)
lines(eps,eps_WR, col="blue",lty=3,lwd=2)
lines(eps,eps_H, col="red",lty=4,lwd=2)
legend(0.5,6,legend=c("M","M+WOR","M+WR","M+H"), col=c("black","green3","blue","red"),
       bty = 'n' ,lwd=c(2,2,2,2),lty=c(1,2,3,4),ncol=1,cex=1.2)

# (b) delta PA vs M plot - Lap
theta=1
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

delta_WR=rep(0,length(eps_WR))
for (i in 1:length(eps_WR)){
  terms=rep(0,n)
  for (k in 1:n){
    delta_k=min(1,(exp(eps[i])-1)*delta_Lap(eps[i]/k)/(exp(eps[i]/k)-1))
    terms[k]=choose(n,k)*(1/n)^k*(1-1/n)^(n-k)*delta_k
  }
  delta_WR[i]=sum(terms)
}

delta_H=rep(0,length(eps_H))
for (i in 1:length(eps_H)){
  terms=rep(0,n)
  for (k in 1:n){
    delta_k=min((exp(eps[i])-1)*delta_Lap(eps[i]/k)/(exp(eps[i]/k)-1),1)
    terms[k]=(b/n)*choose(n,k)*((1/b)^k)*((1-1/b)^(n-k))*delta_k
  }
  delta_H[i]=sum(terms)
}

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

x11()
plot(delta,delta,pch='',xlim=c(0,0.4),ylim=c(0,0.4),xlab=expression(delta),ylab=expression(paste(delta,"'")),cex.lab = 1.5)
lines(delta,delta, col="black",lty=1,lwd=2) # no subsampling
lines(delta,delta_WOR_WB, col="green3",lty=2,lwd=2)
lines(delta,delta_WR_WB, col="blue",lty=3,lwd=2)
lines(delta,delta_H_WB, col="red",lty=4,lwd=2)
legend(0.02,0.4,legend=c("Lap","Lap+WOR","Lap+WR","Lap+H"), col=c("black","green3","blue","red"),
       bty = 'n' ,lwd=c(2,2,2,2),lty=c(1,2,3,4),ncol=1,cex=1.2) 

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

delta_WR=rep(0,length(eps_WR))
for (i in 1:length(eps_WR)){
  terms=rep(0,n)
  for (k in 1:n){
    delta_k=min(1,(exp(eps[i])-1)*delta_Gau(eps[i]/k)/(exp(eps[i]/k)-1))
    terms[k]=choose(n,k)*(1/n)^k*(1-1/n)^(n-k)*delta_k
  }
  delta_WR[i]=sum(terms)
}

delta_H=rep(0,length(eps_H))
for (i in 1:length(eps_H)){
  terms=rep(0,n)
  for (k in 1:n){
    delta_k=min((exp(eps[i])-1)*delta_Gau(eps[i]/k)/(exp(eps[i]/k)-1),1)
    terms[k]=(b/n)*choose(n,k)*((1/b)^k)*((1-1/b)^(n-k))*delta_k
  }
  delta_H[i]=sum(terms)
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
x11()
plot(delta,delta,pch='',xlim=c(0,0.4),ylim=c(0,0.4),xlab=expression(delta),ylab=expression(paste(delta,"'")),cex.lab = 1.5)
lines(delta,delta, col="black",lty=1,lwd=2) # no subsampling
lines(delta,delta_WOR_WB, col="green3",lty=2,lwd=2)
lines(delta,delta_WR_WB, col="blue",lty=3,lwd=2)
lines(delta,delta_H_WB, col="red",lty=4,lwd=2)
legend(0.02,0.4,legend=c("Gauss","Gauss+WOR","Gauss+WR","Gauss+H"), col=c("black","green3","blue","red"),
       bty = 'n' ,lwd=c(2,2,2,2),lty=c(1,2,3,4),ncol=1,cex=1.2) 

#### Figure 2: 6pairs

# (a) Laplace
eps=c(0.05,0.5,1,2,3,4.5)
x11()
plot(delta,eps,pch=4,xlim=c(0,0.5),ylim=c(0,5),ylab=expression(epsilon),xlab=expression(delta),cex.lab = 1.5)
points(delta_WOR_WB,eps_WOR,pch=0)
points(delta_WR_WB,eps_WR,pch=2)
points(delta_H_WB,eps_H,pch=1)
lines(c(delta[1],delta_WOR_WB[1],delta_WR_WB[1],delta_H_WB[1]),c(eps[1],eps_WOR[1],eps_WR[1],eps_H[1]), col="red",lty=1,lwd=2) 
lines(c(delta[2],delta_WOR_WB[2],delta_WR_WB[2],delta_H_WB[2]),c(eps[2],eps_WOR[2],eps_WR[2],eps_H[2]), col="red",lty=1,lwd=2) 
lines(c(delta[3],delta_WOR_WB[3],delta_WR_WB[3],delta_H_WB[3]),c(eps[3],eps_WOR[3],eps_WR[3],eps_H[3]), col="red",lty=1,lwd=2) 
lines(c(delta[4],delta_WOR_WB[4],delta_WR_WB[4],delta_H_WB[4]),c(eps[4],eps_WOR[4],eps_WR[4],eps_H[4]), col="red",lty=1,lwd=2)
lines(c(delta[5],delta_WOR_WB[5],delta_WR_WB[5],delta_H_WB[5]),c(eps[5],eps_WOR[5],eps_WR[5],eps_H[5]), col="red",lty=1,lwd=2) 
lines(c(delta[6],delta_WOR_WB[6],delta_WR_WB[6],delta_H_WB[6]),c(eps[6],eps_WOR[6],eps_WR[6],eps_H[6]), col="red",lty=1,lwd=2) 
legend(0.3,5,legend=c("Lap","Lap+WOR","Lap+WR","Lap+H"), 
       pch = c(4,0,2,1),bty='n')

# (b) Gaussian
eps=c(0.05,0.5,1,2,3,4.5)
x11()
plot(delta,eps,pch=4,xlim=c(0,0.5),ylim=c(0,5),ylab=expression(epsilon),xlab=expression(delta),cex.lab = 1.5)
points(delta_WOR_WB,eps_WOR,pch=0)
points(delta_WR_WB,eps_WR,pch=2)
points(delta_H_WB,eps_H,pch=1)
lines(c(delta[1],delta_WOR_WB[1],delta_WR_WB[1],delta_H_WB[1]),c(eps[1],eps_WOR[1],eps_WR[1],eps_H[1]), col="red",lty=1,lwd=2) 
lines(c(delta[2],delta_WOR_WB[2],delta_WR_WB[2],delta_H_WB[2]),c(eps[2],eps_WOR[2],eps_WR[2],eps_H[2]), col="red",lty=1,lwd=2) 
lines(c(delta[3],delta_WOR_WB[3],delta_WR_WB[3],delta_H_WB[3]),c(eps[3],eps_WOR[3],eps_WR[3],eps_H[3]), col="red",lty=1,lwd=2) 
lines(c(delta[4],delta_WOR_WB[4],delta_WR_WB[4],delta_H_WB[4]),c(eps[4],eps_WOR[4],eps_WR[4],eps_H[4]), col="red",lty=1,lwd=2) 
lines(c(delta[5],delta_WOR_WB[5],delta_WR_WB[5],delta_H_WB[5]),c(eps[5],eps_WOR[5],eps_WR[5],eps_H[5]), col="red",lty=1,lwd=2) 
lines(c(delta[6],delta_WOR_WB[6],delta_WR_WB[6],delta_H_WB[6]),c(eps[6],eps_WOR[6],eps_WR[6],eps_H[6]), col="red",lty=1,lwd=2) 
legend(0.3,5,legend=c("Gauss","Gauss+WOR","Gauss+WR","Gauss+H"), 
       pch = c(4,0,2,1),bty='n')
