n=1000
m=100:150
b=150:200
eta=matrix(0,length(m),length(b))
# OW
for (i in 1:length(m)){
  for (t in 1:length(b)){
    eta[i,t]=(b[t]/n)*(1-(1-1/b[t])^m[i])
  }
}
# WO
eta_WO=matrix(0,length(m),length(b))
for (i in 1:length(m)){
  for (t in 1:length(b)){
    j=seq(1,b[t],1)
    eta_WO[i,t]=sum(choose(b[t],j)*(1/n)^j*(1-1/n)^(b[t]-j)*(1-choose(b[t]-j,m[i])/choose(b[t],m[i])))
  }
}
# WW
eta_WW=matrix(0,length(m),length(b))
for (i in 1:length(m)){
  for (t in 1:length(b)){
    j=seq(1,b[t],1)
    eta_WW[i,t]=sum(choose(b[t],j)*(1/n)^j*(1-1/n)^(b[t]-j)*(1-(1-j/b[t])^m[i]))
  }
}

# eta vs b, m 
# 3D plot
library(plotly)
x11()
plot_ly(x=m,y=b,z=eta, type="surface")%>% layout(scene=list(yaxis=list(title="b"),xaxis =list(title="m"),zaxis=list(title=intToUtf8(0x03B7))))
plot_ly(x=m,y=b,z=eta_WO, type="surface")%>% layout(scene=list(yaxis=list(title="b"),xaxis =list(title="m"),zaxis=list(title=intToUtf8(0x03B7))))
plot_ly(x=m,y=b,z=eta_WW, type="surface")%>% layout(scene=list(yaxis=list(title="b"),xaxis =list(title="m"),zaxis=list(title=intToUtf8(0x03B7))))

# contour plot
plot_ly(x=m,y=b,z=eta, type="contour",contours=list(showlabels=TRUE),line = list(width = 1,color = "white"))%>% layout(yaxis=list(title="b"),xaxis =list(title="m"))
plot_ly(x=m,y=b,z=eta_WO, type="contour",contours=list(showlabels=TRUE),line = list(width = 1,color = "white"))%>% layout(yaxis=list(title="b"),xaxis =list(title="m"))
plot_ly(x=m,y=b,z=eta_WW, type="contour",contours=list(showlabels=TRUE),line = list(width = 1,color = "white"))%>% layout(yaxis=list(title="b"),xaxis =list(title="m"))


# delta'-delta vs b, m, Laplace
theta=1
eps=1
delta_Lap_WB <- function(eps,theta) {
  ifelse(1-exp((eps-theta)/2)<0,0,1-exp((eps-theta)/2))
}
# OW
delta_OW=matrix(0,length(m),length(b))
for (i in 1:length(m)){
  for (t in 1:length(b)){
  terms=rep(0,m[i])
  for (k in 1:m[i]){
    delta_k=delta_Lap_WB(eps,k*theta)
    terms[k]=(b[t]/n)*choose(m[i],k)*((1/b[t])^k)*((1-1/b[t])^(m[i]-k))*delta_k
  }
  delta_OW[i,t]=sum(terms)
}
}
# WO
delta_WO=matrix(0,length(m),length(b))
for (i in 1:length(m)){
  for (t in 1:length(b)){
  terms=rep(0,b[t])
  for (k in 1:m[i]){
    uterms=rep(0,k)
    for (u in 1:k){
      delta_u=delta_Lap_WB(eps,u*theta)
      uterms[u]=(choose(k,u)*choose(b[t]-k,m[i]-u)/choose(b[t],m[i]))*delta_u
    }
    terms[k]=choose(b[t],k)*(1/n)^k*(1-1/n)^(b[t]-k)*sum(uterms)
  }
  for (k in (m[i]+1):b[t]){
    umterms=rep(0,m[i])
    for (u in 1:m[i]){
      delta_u=delta_Lap_WB(eps,u*theta)
      umterms[u]=(choose(k,u)*choose(b[t]-k,m[i]-u)/choose(b[t],m[i]))*delta_u
    }
    terms[k]=choose(b[t],k)*(1/n)^k*(1-1/n)^(b[t]-k)*sum(umterms)
  }
  delta_WO[i,t]=sum(terms)
}
}

delta_WW=matrix(0,length(m),length(b))
for (i in 1:length(m)){
  for (t in 1:length(b)){
  terms=rep(0,b[t])
  for (k in 1:b[t]){
    uterms=rep(0,m[i])
    for (u in 1:m[i]){
      delta_u=delta_Lap_WB(eps,u*theta)
      uterms[u]=choose(m[i],u)*(k/b[t])^u*(1-k/b[t])^(m[i]-u)*delta_u
    }
    terms[k]=choose(b[t],k)*(1/n)^k*(1-1/n)^(b[t]-k)*sum(uterms)
  }
  delta_WW[i,t]=sum(terms)
}
}

x11()
library(plotly)
plot_ly(x=m,y=b,z=delta_OW-delta_Lap_WB(eps,theta), type="surface")%>% layout(scene=list(yaxis=list(title="b"),xaxis =list(title="m"),zaxis=list(title=paste(intToUtf8(0x03B4),"'-",intToUtf8(0x03B4)))))
plot_ly(x=m,y=b,z=delta_WO-delta_Lap_WB(eps,theta), type="surface")%>% layout(scene=list(yaxis=list(title="b"),xaxis =list(title="m"),zaxis=list(title=paste(intToUtf8(0x03B4),"'-",intToUtf8(0x03B4)))))
plot_ly(x=m,y=b,z=delta_WW-delta_Lap_WB(eps,theta), type="surface")%>% layout(scene=list(yaxis=list(title="b"),xaxis =list(title="m"),zaxis=list(title=paste(intToUtf8(0x03B4),"'-",intToUtf8(0x03B4)))))

# contour plot
plot_ly(x=m,y=b,z=delta_OW-delta_Lap_WB(eps,theta), type="contour",contours=list(showlabels=TRUE),line = list(width = 1,color = "white"))%>% layout(yaxis=list(title="b"),xaxis =list(title="m"))
plot_ly(x=m,y=b,z=delta_WO-delta_Lap_WB(eps,theta), type="contour",contours=list(showlabels=TRUE),line = list(width = 1,color = "white"))%>% layout(yaxis=list(title="b"),xaxis =list(title="m"))
plot_ly(x=m,y=b,z=delta_WW-delta_Lap_WB(eps,theta), type="contour",contours=list(showlabels=TRUE),line = list(width = 1,color = "white"))%>% layout(yaxis=list(title="b"),xaxis =list(title="m"))



# delta'-delta vs b, m, Gaussian
theta=1
eps=1
delta_Gau_WB <- function(eps,theta) {
  del=pnorm(theta/2-eps/theta)-exp(eps)*pnorm(-theta/2-eps/theta)
  return (min(del,1))
}

# OW
delta_OW_Gauss=matrix(0,length(m),length(b))
for (i in 1:length(m)){
  for (t in 1:length(b)){
    terms=rep(0,m[i])
    for (k in 1:m[i]){
      delta_k=delta_Gau_WB(eps,k*theta)
      terms[k]=(b[t]/n)*choose(m[i],k)*((1/b[t])^k)*((1-1/b[t])^(m[i]-k))*delta_k
    }
    delta_OW_Gauss[i,t]=sum(terms)
  }
}
# WO
delta_WO_Gauss=matrix(0,length(m),length(b))
for (i in 1:length(m)){
  for (t in 1:length(b)){
    terms=rep(0,b[t])
    for (k in 1:m[i]){
      uterms=rep(0,k)
      for (u in 1:k){
        delta_u=delta_Gau_WB(eps,u*theta)
        uterms[u]=(choose(k,u)*choose(b[t]-k,m[i]-u)/choose(b[t],m[i]))*delta_u
      }
      terms[k]=choose(b[t],k)*(1/n)^k*(1-1/n)^(b[t]-k)*sum(uterms)
    }
    for (k in (m[i]+1):b[t]){
      umterms=rep(0,m[i])
      for (u in 1:m[i]){
        delta_u=delta_Gau_WB(eps,u*theta)
        umterms[u]=(choose(k,u)*choose(b[t]-k,m[i]-u)/choose(b[t],m[i]))*delta_u
      }
      terms[k]=choose(b[t],k)*(1/n)^k*(1-1/n)^(b[t]-k)*sum(umterms)
    }
    delta_WO_Gauss[i,t]=sum(terms)
  }
}
# WW
delta_WW_Gauss=matrix(0,length(m),length(b))
for (i in 1:length(m)){
  for (t in 1:length(b)){
    terms=rep(0,b[t])
    for (k in 1:b[t]){
      uterms=rep(0,m[i])
      for (u in 1:m[i]){
        delta_u=delta_Gau_WB(eps,u*theta)
        uterms[u]=choose(m[i],u)*(k/b[t])^u*(1-k/b[t])^(m[i]-u)*delta_u
      }
      terms[k]=choose(b[t],k)*(1/n)^k*(1-1/n)^(b[t]-k)*sum(uterms)
    }
    delta_WW_Gauss[i,t]=sum(terms)
  }
}

x11()
library(plotly)
plot_ly(x=m,y=b,z=delta_OW_Gauss-delta_Gau_WB(eps,theta), type="surface")%>% layout(scene=list(yaxis=list(title="b"),xaxis =list(title="m"),zaxis=list(title=paste(intToUtf8(0x03B4),"'-",intToUtf8(0x03B4)))))
plot_ly(x=m,y=b,z=delta_WO_Gauss-delta_Gau_WB(eps,theta), type="surface")%>% layout(scene=list(yaxis=list(title="b"),xaxis =list(title="m"),zaxis=list(title=paste(intToUtf8(0x03B4),"'-",intToUtf8(0x03B4)))))
plot_ly(x=m,y=b,z=delta_WW_Gauss-delta_Gau_WB(eps,theta), type="surface")%>% layout(scene=list(yaxis=list(title="b"),xaxis =list(title="m"),zaxis=list(title=paste(intToUtf8(0x03B4),"'-",intToUtf8(0x03B4)))))

# contour plot
plot_ly(x=m,y=b,z=delta_OW_Gauss-delta_Gau_WB(eps,theta), type="contour",contours=list(showlabels=TRUE),line = list(width = 1,color = "white"))%>% layout(yaxis=list(title="b"),xaxis =list(title="m"))
plot_ly(x=m,y=b,z=delta_WO_Gauss-delta_Gau_WB(eps,theta), type="contour",contours=list(showlabels=TRUE),line = list(width = 1,color = "white"))%>% layout(yaxis=list(title="b"),xaxis =list(title="m"))
plot_ly(x=m,y=b,z=delta_WW_Gauss-delta_Gau_WB(eps,theta), type="contour",contours=list(showlabels=TRUE),line = list(width = 1,color = "white"))%>% layout(yaxis=list(title="b"),xaxis =list(title="m"))


