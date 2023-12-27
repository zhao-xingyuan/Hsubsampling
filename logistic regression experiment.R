rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load the preprocessed data
adult_original=read.csv("C:/Users/Zhao_/Desktop/RUB/bootstrap/adult_processed2.csv",header=TRUE)
adult=adult_original
for(j in 1:ncol(adult_original)){
  adult[,j]=adult_original[,j]/max(adult_original[,j])
}
test_size=nrow(adult)/4  #10323
train_size=nrow(adult)*(3/4) #30969

adult_train=adult[1:train_size,]
adult_test=adult[(train_size+1):nrow(adult),]

y_train=adult_train$income
table(y_train)
#0    1 
#23203  7766 
X_train=adult_train[,-23]
y_test=adult_test$income
table(y_test)
#0   1 
#7641 2682  
X_test=adult_test[,-23]
test_matrix = data.matrix(cbind(Intercept = 1, X_test))

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

sigmoid<- function(z){
  1/(1+exp(-z))
}

sgd <- function(
  stepsize = 0.0008, regularization = FALSE, lambda=0.1,
  subsample ='full-non private', C=3, TT= 200, eps_prime=0.1, n=train_size, b= 500, m= 300)
{
  n  = nrow(adult_train)
  y  = y_train
  X  = data.matrix(cbind(Intercept = 1, X_train))
  train=cbind(X,y)
  
  compute_sigma <- function(sensitivity, eps){
    sqrt(2*log(1.25/(1/n)))*sensitivity/(eps*n)
  }
  
  if(subsample=='rub') {
    eps_H=log(1+(exp(eps_prime)-1)/((b/n)*(1-(1-1/b)^m)))
    sigma=compute_sigma(C, eps_H)
  }
  if(subsample=='full-non private') {
    sigma=0
  }
  if(subsample=='WO') {
    WO_terms=rep(0,b)
    for (j in 1:b){
      WO_terms[j]=choose(b,j)*(1/n)^j*(1-1/n)^(b-j)*(1-choose(b-j,m)/choose(b,m))
    }
    eta_WO=sum(WO_terms)
    eps_WO=log(1+(exp(eps_prime)-1)/eta_WO)
    sigma=compute_sigma(C, eps_WO)
  }
  if(subsample=='WW') {
    WW_terms=rep(0,b)
    for (j in 1:b){
      WW_terms[j]=choose(b,j)*(1/n)^j*(1-1/n)^(b-j)*(1-(1-j/b)^m)
    }
    eta_WW=sum(WW_terms)
    eps_WW=log(1+(exp(eps_prime)-1)/eta_WW)
    sigma=compute_sigma(C, eps_WW)
  }
  if(subsample=='WR') {
    eps_WR=log(1+(exp(eps_prime)-1)/(1-(1-1/n)^m))
    sigma=compute_sigma(C, eps_WR)
  }
  if(subsample=='WOR') {
    eps_WOR=log(1+(exp(eps_prime)-1)/(m/n))
    sigma=compute_sigma(C, eps_WOR)
  }
  if(subsample=='Poisson') {
    sigma=0
  }
  if(subsample=='DP-Poisson') {
    eps_Poisson=log(1+(exp(eps_prime)-1)/(m/n))
    sigma=compute_sigma(C, eps_Poisson)
  }
  
  bb = rep(1,ncol(X))  # initial value of the parameters
  loss = rep(0,TT)                                    
  for (t in 1:TT) {
    if(subsample=='rub') {
      s= RUB(1:n,b,m)
    }
    if(subsample=='full-non private') {
      s= seq(1,n,1) 
    }
    if(subsample=='WO') {
      s= WO(1:n,b,m)
    }
    if(subsample=='WW') {
      s= RUB(1:n,b,m)
    }
    if(subsample=='WR') {
      s= WR(1:n,m)
    }
    if(subsample=='WOR') {
      s= WOR(1:n,m)  
    }
    if(subsample=='Poisson') {
      s= Poisson(1:n,m/n) 
    }
    if(subsample=='DP-Poisson') {
      s= Poisson(1:n,m/n) 
    }
    grad= matrix(0,length(s),ncol(X))
    
    for(i in 1:length(s)){
      Xi   = X[s[i], ]; 
      yi   = y[s[i]]; 
      if (regularization == FALSE){
        gradi = -Xi%*%(yi*(exp(Xi%*%bb)+1)^(-1))+Xi%*%((1-yi)*(exp(-Xi%*%bb)+1)^(-1))
      }else{
        gradi = -Xi%*%(yi*(exp(Xi%*%bb)+1)^(-1))+Xi%*%((1-yi)*(exp(-Xi%*%bb)+1)^(-1))+2*lambda*bb
      }

      if (subsample=='Poisson'|subsample=='full-non private') {
        grad[i,] = t(gradi)
      }
      else{grad[i,] = t(gradi/max(1,norm(gradi, type="2")/C))}
    }
    #grads=apply(grad,2,sum)/length(s)+rnorm(ncol(X),0, sigma)
    grads=(apply(grad,2,sum)+(length(s)/m)*rnorm(ncol(X),0, m*sigma))/length(s)
    bb = bb - stepsize*grads
    if (regularization == FALSE){
      loss[t] = -mean(y*as.vector(log(sigmoid(X%*%bb)))+(1-y)*as.vector(log(1-sigmoid(X%*%bb))))
    }else{
      loss[t] = -mean(y*as.vector(log(sigmoid(X%*%bb)))+(1-y)*as.vector(log(1-sigmoid(X%*%bb))))+lambda*(norm(bb, type="2"))^2
    }
  }
  # prediction
  log_odds=test_matrix%*%bb
  predict_y=ifelse(log_odds>0,1,0)
  if (sum(predict_y)==0){
    confusion=rbind(table(predict_y,y_test),c(0,0))
  }
  else if(sum(predict_y)==length(predict_y)){
    confusion=rbind(c(0,0),table(predict_y,y_test))
  }
  else{confusion=table(predict_y,y_test)}
  metrics=rep(0,6) #accuracy	sensitivity	specificity	PPV 	NPV	F1-score
  metrics[1]=(confusion[1,1]+confusion[2,2])/length(y_test)
  metrics[2]=confusion[2,2]/(confusion[1,2]+confusion[2,2])
  metrics[3]=confusion[1,1]/(confusion[1,1]+confusion[2,1])
  metrics[4]=confusion[2,2]/(confusion[2,1]+confusion[2,2])
  metrics[5]=confusion[1,1]/(confusion[1,1]+confusion[1,2])
  metrics[6]=2*metrics[2]*metrics[4]/(metrics[2]+metrics[4])

  list(est= bb,loss=loss, sig=sigma, cross_entropy=mean(loss[(TT-15):TT]),prediction=metrics)
}


# eps'=0.0005, T=1000, lr=0.01
metric=matrix(0,7,6)
TT=1000
loss_history=matrix(0,7,TT)
noise_sig=rep(0,7)
lr=0.01
eps_p=0.0005
set.seed(9)
x11()
start_time<-Sys.time()
out=sgd(subsample="WR", TT=TT,stepsize=lr, C=3,eps_prime=eps_p)
metric[2,]=out$prediction
loss_history[2,]=out$loss
noise_sig[2]=out$sig
plot(1:TT,loss_history[2,],type="l", xlab="", ylab="", ylim=c(0,3.5), col='green3',lty = 'aa', lwd = 2) 
out=sgd(subsample="WOR", TT=TT, stepsize=lr, C=3, eps_prime=eps_p) 
metric[3,]=out$prediction 
loss_history[3,]=out$loss
noise_sig[3]=out$sig
lines(1:TT,loss_history[3,],type="l", col='blue',lty = 2, lwd = 2)
out=sgd(subsample="rub", TT=TT, C=3, stepsize=lr, eps_prime=eps_p)
metric[5,]=out$prediction
loss_history[5,]=out$loss
noise_sig[5]=out$sig
lines(1:TT,loss_history[5,],type="l", col='red',lty = 3, lwd = 2)
out=sgd(subsample="WO", TT=TT, C=3, stepsize=lr, eps_prime=eps_p)
metric[6,]=out$prediction
loss_history[6,]=out$loss
noise_sig[6]=out$sig
lines(1:TT,loss_history[6,],type="l", col='orange',lty = 4, lwd = 2)
out=sgd(subsample="WW", TT=TT, C=3, stepsize=lr, eps_prime=eps_p)
metric[7,]=out$prediction
loss_history[7,]=out$loss
noise_sig[7]=out$sig
lines(1:TT,loss_history[7,],type="l", col='purple',lty = 5, lwd = 2)
out=sgd(subsample="DP-Poisson", TT=TT, C=3, stepsize=lr, eps_prime=eps_p)
metric[4,]=out$prediction
loss_history[4,]=out$loss
noise_sig[4]=out$sig
lines(1:TT,loss_history[4,],type="l", col='brown',lty = 6, lwd = 2)
out_nonprivate=sgd(subsample="Poisson", TT=TT, C=3,stepsize=lr,  eps_prime=eps_p)
metric[1,]=out_nonprivate$prediction 
loss_history[1,]=out_nonprivate$loss
noise_sig[1]=out$sig
lines(1:TT,loss_history[1,],type="l", col='black',lty = 1, lwd = 2)
legend(300, 3.1, c("Poisson (non-private)","Poisson","WOR", "WR", "MUST.OW","MUST.WW"),  col=c('black','brown','green3','blue','red','purple'),lty=c('solid','twodash','aa','dashed','dotted','longdash'),bty = "n")
mtext("loss",2,2)
mtext("iteration",1,2,cex=1)
end_time<-Sys.time()
end_time-start_time #3.62864 mins

m=300
write.csv(noise_sig*m, "noise sigma_eps0005.csv")
write.csv(metric, "utility logistic metrics_T1000_eps0005.csv")
write.csv(loss_history, "loss_history_T1000_eps0005.csv")
