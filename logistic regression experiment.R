###########################################################################
#### Logistic Regression DP-SGD 
###########################################################################



###########################################################################
#### 1. Preparation
###########################################################################
# load the preprocessed data
adult_original=read.csv("adult_processed2.csv",header=TRUE)
adult=adult_original
for(j in 1:ncol(adult_original)){
  adult[,j]=adult_original[,j]/max(adult_original[,j])
}


# Split train and test
test_size=floor(nrow(adult)/4)  #10323
train_size=floor(nrow(adult)*(3/4)) #30969
adult_train=adult[1:train_size,]
adult_test=adult[(train_size+1):(train_size+test_size),]

y_train=adult_train$income
table(y_train)
#0    1 
#23203  7766 
X_train=adult_train[,-14]
y_test=adult_test$income
table(y_test)
#0   1 
#7641 2682  
X_test=adult_test[,-14]
test_matrix = data.matrix(cbind(Intercept = 1, X_test))
n=train_size



###########################################################################
#### 2. Functions
###########################################################################

#-------------------------------------------------------------------
# 2.1 Functions for computing sigma for Gaussian mechanism in DP-SGD 
#-------------------------------------------------------------------
compute_sigma_original <- function(sensitivity, eps){
  sqrt(2*log(1.25/(1/n)))*sensitivity/(eps*n)
}
compute_sigma <- function(sensitivity, eps){
  delta <- 1/n
  delta_function <- function(sigma) {
    theta <- sensitivity / (n * sigma)
    pnorm(-eps/theta + theta/2) - exp(eps) * pnorm(-eps/theta - theta/2) - delta
  }
  sigma_root = uniroot(delta_function, c(1e-10, 100))$root
  return(sigma_root)
}

# For check, compute_sigma is tighter
for(i in seq(0.01, 1, 0.05)){
  print(c(compute_sigma_original(3, i), compute_sigma(3,i)))
}

#-------------------------------------------------------------------
# 2.2 Functions for different sampling schemes 
#-------------------------------------------------------------------
OW = function(x, b, m){
  xb = sample(x, b, replace=F)
  xOW = sample(xb,m, replace=T)
  return(xOW)
}

WO = function(x, b, m){
  xb = sample(x, b, replace=T)
  xWO = sample(xb,m, replace=F)
  return(xWO)
}

WW = function(x, b, m){
  xb = sample(x, b, replace=T)
  xWW = sample(xb,m, replace=T)
  return(xWW)
}

WR  =  function(x, m)  sample(x, m, replace=T)
WOR =  function(x, m)  sample(x, m, replace=F)

Poisson = function(x,q){
  random_values <- runif(length(x))
  xpoisson <- x[random_values < q]
  return(xpoisson)
}

sigmoid<- function(z){
  1/(1+exp(-z))
}


#-------------------------------------------------------------------
# 2.3 Function for DP-SGD
#-------------------------------------------------------------------
sgd <- function(
    stepsize = 0.0005, regularization = FALSE, lambda=0.1,
    subsample ='full-non private', C=3, TT= 200, eps_prime=0.1, n=train_size, b= 10, m= 10)
{
  n  = train_size
  y  = y_train
  X  = data.matrix(cbind(Intercept = 1, X_train))
  train=cbind(X,y)

  # Calculate per-query eps and corresponding delta
  if(subsample=='OW') {
    eps_H=log(1+(exp(eps_prime)-1)/((b/n)*(1-(1-1/b)^m)))
    sigma=compute_sigma(C, eps_H)
  }else if(subsample=='WO') {
    j_values <- 1:b  
    
    choose_b_j <- choose(b, j_values)  # choose(b, j)
    term1 <- (1/n)^j_values            # (1/n)^j
    term2 <- (1 - 1/n)^(b - j_values)  # (1 - 1/n)^(b - j)
    
    WO_terms <- choose_b_j * term1 * term2 * (1 - choose(b - j_values, m) / choose(b, m))
    eps_WO=log(1+(exp(eps_prime)-1)/sum(WO_terms))
    sigma=compute_sigma(C, eps_WO)
  }else if(subsample=='WW') {
    j_values <- 1:b  
    
    choose_b_j <- choose(b, j_values)  # choose(b, j)
    term1 <- (1/n)^j_values            # (1/n)^j
    term2 <- (1 - 1/n)^(b - j_values)  # (1 - 1/n)^(b - j)
    term3 <- (1 - (1 - j_values / b)^m) # (1 - (1 - j/b)^m)
    
    WW_terms <- choose_b_j * term1 * term2 * term3
    eps_WW=log(1+(exp(eps_prime)-1)/sum(WW_terms))
    sigma=compute_sigma(C, eps_WW)
  }else if(subsample=='WR') {
    eps_WR=log(1+(exp(eps_prime)-1)/(1-(1-1/n)^m))
    sigma=compute_sigma(C, eps_WR)
  }else if(subsample=='WOR') {
    eps_WOR=log(1+(exp(eps_prime)-1)/(m/n))
    sigma=compute_sigma(C, eps_WOR)
  }else if(subsample=='Poisson') {
    sigma=0
  }else if(subsample=='full-non private'){
    sigma=0
  }else if(subsample=='DP-Poisson') {
    eps_Poisson=log(1+(exp(eps_prime)-1)/(m/n))
    sigma=compute_sigma(C, eps_Poisson)
  }
  
  bb = rep(1,ncol(X))  # initial value of the parameters
  loss = rep(0,TT) 
  grad_diff = rep(0, TT)
  for (t in 1:TT) {
    # Subsampling
    if(subsample=='OW') {
      s= OW(1:n,b,m)
    }else if(subsample=='WO') {
      s= WO(1:n,b,m)
    }else if(subsample=='WW') {
      s= WW(1:n,b,m)
    }else if(subsample=='WR') {
      s= WR(1:n,m)
    }else if(subsample=='WOR') {
      s= WOR(1:n,m)  
    }else if(subsample=='Poisson') {
      s= Poisson(1:n,m/n) 
    }
    else if(subsample=='full-non private') {
      s= seq(1,n,1) 
    }else if(subsample=='DP-Poisson') {
      s= Poisson(1:n,m/n) 
    }
    
    # Gradient Descent
    unique_s = unique(s)
    freq_s = table(s)
    grad= matrix(0,length(unique_s),ncol(X))
    grad_ori= matrix(0,length(unique_s),ncol(X))
    
    # Compute the gradient only for the unique element
    for(i in seq_along(unique_s)){
      idx <- unique_s[i]

      Xi = X[idx,]
      yi <- y[idx]  
      if (regularization == FALSE){
        gradi = -Xi%*%(yi*(exp(Xi%*%bb)+1)^(-1))+Xi%*%((1-yi)*(exp(-Xi%*%bb)+1)^(-1))
      }else{
        gradi = -Xi%*%(yi*(exp(Xi%*%bb)+1)^(-1))+Xi%*%((1-yi)*(exp(-Xi%*%bb)+1)^(-1))+2*lambda*bb
      }
      
      # Apply gradient clipping and noise if necessary
      if (subsample == 'Poisson' || subsample == 'full-non private') {
        grad[i, ] <- t(gradi)
      } else{
        grad_ori[i, ] <- t(gradi / max(1, norm(gradi, type = "2") / C))
        grad[i, ] <- grad_ori[i, ] + rnorm(ncol(X), 0, sigma)
      }
    }
    grads <- colSums(grad * as.numeric(freq_s)) / length(s)
    grads_ori <- colSums(grad_ori * as.numeric(freq_s)) / length(s)
    grad_diff[t] = grads[1] - grads_ori[1]
    
    # print(c(grads_diff = grads[1] - grads_ori[1], bb_diff = (bb - stepsize*grads)[1]-(bb - stepsize*grads_ori)[1]))    
    
    
    # Update estimate
    bb = bb - stepsize*grads

    
    # Calculate loss
    if (regularization == FALSE){
      loss[t] = -mean(y*as.vector(log(sigmoid(X%*%bb)))+(1-y)*as.vector(log(1-sigmoid(X%*%bb))), na.rm = TRUE)
      # loss[t] = mean((y - sigmoid(X %*% bb))^2, na.rm = TRUE)
    }else{
      loss[t] = -mean(y*as.vector(log(sigmoid(X%*%bb)))+(1-y)*as.vector(log(1-sigmoid(X%*%bb))))+lambda*(norm(bb, type="2"))^2
    }
  }
  
  # prediction
  log_odds=test_matrix%*%bb
  predict_y=ifelse(log_odds>0,1,0)
  if (sum(predict_y)==0){
    confusion=rbind(table(predict_y,y_test),c(0,0))
  }else if(sum(predict_y)==length(predict_y)){
    confusion=rbind(c(0,0),table(predict_y,y_test))
  }else{confusion=table(predict_y,y_test)}
  
  metrics = rep(0,6) #accuracy	sensitivity	specificity	PPV 	NPV	F1-score
  metrics[1]=(confusion[1,1]+confusion[2,2])/length(y_test)
  metrics[2]=confusion[2,2]/(confusion[1,2]+confusion[2,2])
  metrics[3]=confusion[1,1]/(confusion[1,1]+confusion[2,1])
  metrics[4]=confusion[2,2]/(confusion[2,1]+confusion[2,2])
  metrics[5]=confusion[1,1]/(confusion[1,1]+confusion[1,2])
  metrics[6]=2*metrics[2]*metrics[4]/(metrics[2]+metrics[4])
  names(metrics) <- c("Accuracy", "Sensitivity", "Specificity", "PPV", "NPV", "F1-score")
  
  list(est= bb,loss=loss, sig=sigma, grad_diff = grad_diff,
       cross_entropy=mean(loss[1:TT]),prediction=metrics)
}


combinations = expand.grid(eps=c(5e-10, 5e-4, 5e-1, 5), C=c(0.5, 1, 2, 3))
for(i in 1:nrow(combinations)){
  set.seed(1)
  combinations$avg_grad_diff_b500_m300[i]=formatC(mean(sgd(subsample="WR", TT=20,stepsize=1, C=combinations$C[i],eps_prime=combinations$eps[i], b=500, m=300)$grad_diff), format = "e", digits = 0)
  combinations$avg_grad_diff_b500_m100[i]=formatC(mean(sgd(subsample="WR", TT=20,stepsize=1, C=combinations$C[i],eps_prime=combinations$eps[i], b=500, m=100)$grad_diff), format = "e", digits = 0)
  combinations$avg_grad_diff_b500_m30[i]=formatC(mean(sgd(subsample="WR", TT=20,stepsize=1, C=combinations$C[i],eps_prime=combinations$eps[i], b=500, m=30)$grad_diff), format = "e", digits = 0)
  
  combinations$avg_grad_diff_b100_m100[i]=formatC(mean(sgd(subsample="WR", TT=20,stepsize=1, C=combinations$C[i],eps_prime=combinations$eps[i], b=100, m=100)$grad_diff), format = "e", digits = 0)
  combinations$avg_grad_diff_b100_m30[i]=formatC(mean(sgd(subsample="WR", TT=20,stepsize=1, C=combinations$C[i],eps_prime=combinations$eps[i], b=100, m=30)$grad_diff), format = "e", digits = 0)
  combinations$avg_grad_diff_b50_m30[i]=formatC(mean(sgd(subsample="WR", TT=20,stepsize=1, C=combinations$C[i],eps_prime=combinations$eps[i], b=50, m=30)$grad_diff), format = "e", digits = 0)
}






sgd(subsample="WW", TT=20,stepsize=0.1, C=2,eps_prime=1e-10, b= b, m=m)
sgd(subsample="Poisson", TT=20,stepsize=1, C=1,eps_prime=1e-10, b= b, m=m)

###########################################################################
#### 3. Plot
###########################################################################

C=1.5
eps_prime = 1e-5
b=200
m=100
TT=15
lr=0.4

metric=as.data.frame(matrix(0,7,6))
loss_history=as.data.frame(matrix(0,7,TT))
noise_sig=rep(0,7)
time_record = rep(0,7)

set.seed(48)
start_time<-Sys.time()
out=sgd(subsample="WR", TT=TT,stepsize=lr, C=C,eps_prime=eps_p, b= b, m=m)
end_time<-Sys.time()
loss_history[2,]=out$loss
metric[2,]=out$prediction
noise_sig[2]=out$sig
time_record[2] = end_time - start_time
par(cex.axis=1.4) 
plot(1:TT,loss_history[2,],type="l", xlab="", ylab="", ylim=c(0,4), col='green3',lty = 'aa', lwd = 2) 

start_time<-Sys.time()
out_nonprivate=sgd(subsample="Poisson", TT=TT, C=C,stepsize=lr,  eps_prime=eps_p, b= b, m=m)
end_time<-Sys.time()
metric[1,]=out_nonprivate$prediction 
loss_history[1,]=out_nonprivate$loss
noise_sig[1]=out_nonprivate$sig
time_record[1] = end_time - start_time
lines(1:TT,loss_history[1,],type="l", col='black',lty = 1, lwd = 2)

start_time<-Sys.time()
out=sgd(subsample="WOR", TT=TT, stepsize=lr, C=C, eps_prime=eps_p, b= b, m=m)
end_time<-Sys.time()
metric[3,]=out$prediction 
loss_history[3,]=out$loss
noise_sig[3]=out$sig
time_record[3] = end_time - start_time
lines(1:TT,loss_history[3,],type="l", col='blue',lty = 2, lwd = 2)

start_time<-Sys.time()
out=sgd(subsample="DP-Poisson", TT=TT, C=C, stepsize=lr, eps_prime=eps_p, b= b, m=m)
end_time<-Sys.time()
metric[4,]=out$prediction
loss_history[4,]=out$loss
noise_sig[4]=out$sig
time_record[4] = end_time - start_time
lines(1:TT,loss_history[4,],type="l", col='brown',lty = 6, lwd = 2)

start_time<-Sys.time()
out=sgd(subsample="OW", TT=TT, C=C, stepsize=lr, eps_prime=eps_p, b= b, m=m)
end_time<-Sys.time()
metric[5,]=out$prediction
loss_history[5,]=out$loss
noise_sig[5]=out$sig
time_record[5] = end_time - start_time
lines(1:TT,loss_history[5,],type="l", col='red',lty = 3, lwd = 2)

# out=sgd(subsample="WO", TT=TT, C=C, stepsize=lr, eps_prime=eps_p, b= b, m=m)
# metric[6,]=out$prediction
# loss_history[6,]=out$loss
# noise_sig[6]=out$sig
# lines(1:TT,loss_history[6,],type="l", col='orange',lty = 4, lwd = 2)

start_time<-Sys.time()
out=sgd(subsample="WW", TT=TT, C=C, stepsize=lr, eps_prime=eps_p, b= b, m=m)
end_time<-Sys.time()
metric[7,]=out$prediction
loss_history[7,]=out$loss
noise_sig[7]=out$sig
time_record[7] = end_time - start_time
lines(1:TT,loss_history[7,],type="l", col='purple',lty = 5, lwd = 2)



rownames(metric) = rownames(loss_history)= c("Poisson", "WR", "WOR", "DP-Poisson", "OW", "WO", "WW")
names(noise_sig) = names(time_record) = c("Poisson", "WR", "WOR", "DP-Poisson", "OW", "WO", "WW")

legend(8, 3.5, c("Poisson (non-private)","Poisson","WOR", "WR", "MUST.ow","MUST.ww"), lwd = 2, seg.len = 1.4, cex = 1.3, 
       col=c('black','brown','green3','blue','red','purple'),lty=c('solid','twodash','aa','dashed','dotted','longdash'),bty = "n")
mtext("loss",2,2, cex = 1.7)
mtext("iteration",1,2,cex=1.7)

write.csv(time_record, "time_record.csv")
write.csv(loss_history, "loss_history.csv")
write.csv(noise_sig, "noise_sig.csv")
colnames(metric) <- c("Accuracy", "Sensitivity", "Specificity", "PPV", "NPV", "F1-score")
write.csv(metric, "metric.csv")
print(metric)

