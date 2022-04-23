rm(list=ls())

beta0 = 0.5
N=100  
library(survival)


#par(mfrow=c(1,2))

condional_haz = function(Covar,times,del,t,x){
  
  n= length(times)
  K = function(x) {1/sqrt(2*pi) * exp(-(x^2)/2)}
  W = function(y) { integrate(K,lower=-Inf,upper=y)$value  }
  
  #The optimal bandwidths are both of order n^(???1/6)
  #and the corresponding AMISE is of order n^(???2/3)
  
  wx=rep(0,n)
  sumw = 0
  for (i in 1:n){
    sumw = sumw+K((x-Covar[i])/hx)
  }
  sumdeno = 0
  sumnumer = 0
  
  for (i in 1:n){
    wx[i] = K((x-Covar[i])/hx)/sumw
    sumdeno = sumdeno+wx[i]*W((times[i]-t)/ht)
    sumnumer = sumnumer+wx[i]*K((t-times[i])/ht)*del[i]
  }
  
  sumnumer/sumdeno/ht
  
}

kernelreg<-function(x,y,t){
  n <- length(x)
  K <- function(x)  1/sqrt(2*pi) * exp(-(x^2)/2)
  h= n^(-0.2)  #Silverman's Rule of Thumb
  numer=0
  deno = 0
  for(i in 1:n){
    numer = numer+y[i]*K((t-x[i])/h)
    deno= deno+K((t-x[i])/h)
  }
  return (numer/deno)
  
}



inversion_gen<-function(n,Z){
  temp = runif(n)
  x = (-1+sqrt(1-40*log(1-temp)/exp(beta0*Z)))/2
}

MHD1 = function(b){
  zeroset= rep(0,length(alltime))
  est = rep(0,length(alltime))
  for (i in 1:length(alltime)){
    zeroset[i]=condional_haz(Z,alltime,del,alltime[i],Z[i])/exp(b*Z[i])
    est[i] = condional_haz(Z,alltime,del,alltime[i],0)
  }
  cons_h0 = function(t){kernelreg(alltime,zeroset,t)}
  #  est_h0= function(t){condional_haz(Z,alltime,del,t,0)}
  # est_h0 = function(t) 0.1+0.2*t
  distance= 0
  for (i in 1:length(alltime)){
    distance = distance+(cons_h0(alltime[i])^0.5-est[i]^0.5)^2
  }
  distance
}

MHD2 = function(b){
  zeroset= rep(0,length(alltime))
  est = rep(0,length(alltime))
  for (i in 1:length(alltime)){
    zeroset[i]=condional_haz(Z,alltime,del,alltime[i],Z[i])/exp(b*Z[i])
    est[i] = condional_haz(Z,alltime,del,alltime[i],0)
  }
  #  cons_h0 = function(t){kernelreg(alltime,zeroset,t)}
  #   est_h0= function(t){condional_haz(Z,alltime,del,t,0)}
  # est_h0 = function(t) 0.1+0.2*t
  distance= 0
  for (i in 1:length(alltime)){
    distance = distance+(zeroset[i]^0.5-est[i]^0.5)^2
  }
  distance
}


MHD3 = function(b){
  zeroset= rep(0,length(alltime))
  # est = rep(0,length(alltime))
  for (i in 1:length(alltime)){
    zeroset[i]=condional_haz(Z,alltime,del,alltime[i],Z[i])/exp(b*Z[i])
    #  est[i] = condional_haz(Z,alltime,del,alltime[i],1)/exp(b)
  }
  cons_h0 = function(t){kernelreg(alltime,zeroset,t)}
  # est_h0= function(t){condional_haz(Z,alltime,del,t,1)}
  est_h0 = function(t) 0.1+0.2*t
  distance= 0
  for (i in 1:length(alltime)){
    distance = distance+(zeroset[i]^0.5-est_h0(alltime[i])^0.5)^2
  }
  distance
}
#est_h0= function(t){condional_haz(Z,alltime,del,t,0)}
#est0 = rep(0,N)
#for (i in 1:N){
#  est0[i] = est_h0(alltime[i])
#}
#zeroset= rep(0,length(alltime))
#for (i in 1:length(alltime)){
#  zeroset[i]=condional_haz(Z,alltime,del,alltime[i],Z[i])/exp(beta0*Z[i])
#}
#bs = function(t) 0.1+0.2*t
#par(mfrow=c(1,3))

#plot(alltime,bs(alltime))
#plot(alltime,est0)
#plot(alltime,zeroset)
#const = function(t) kernelreg(alltime,zeroset,t)
#plot(alltime,const(alltime))
#dis = function(a,b){
#  dis = 0
#  for (i in length(a)){
#    dis = dis+(a[i]^0.5-b[i]^0.5)^2
#  }
#  dis
#}
#dis(zeroset,bs(alltime))
#dis(zeroset,const(alltime))


outcomes1 = c()
outcomes2 = c()
outcomes3 = c()
outcomes4 = c()
outcomes5 = c()
outcomes6 = c()
ptm<-proc.time()
for (i in 1:100){
  Z = runif(N,0,2)
  alltime= c()
  for (i in 1:N){
    temp = inversion_gen(1,Z[i])
    alltime = c(alltime,temp)
  }
  del=rep(1,N)
  # 20% censor
  for (i in 1:N){
    temp = runif(1,0,9)
    if (temp<alltime[i]){
      alltime[i]=temp
      del[i]=0
    }
  }

  
  hx=N^(-1/5)*sd(Z)
  ht=N^(-1/5)*sd(alltime)
  hat=optimize(MHD2,lower = 0,upper = 10)
  outcomes1 = c(outcomes1 , hat$minimum)
  
  hat=optimize(MHD3,lower = 0,upper = 10)
  outcomes2 = c(outcomes2 , hat$minimum)
  
  fit = coxph(Surv(alltime,del)~Z)
  outcomes3 = c(outcomes3,fit$coefficients[1])
  
  #################################################
  index = sample(which(del==1),1)
  alltime[index] = 4
  Z[index] = 4
  
  
  hat=optimize(MHD2,lower = 0,upper = 10)
  outcomes4 = c(outcomes4 , hat$minimum)
  
  hat=optimize(MHD3,lower = 0,upper = 10)
  outcomes5 = c(outcomes5 , hat$minimum)
  
  fit = coxph(Surv(alltime,del)~Z)
  outcomes6 = c(outcomes6,fit$coefficients[1])
} 
mean(outcomes2)- beta0   # MHDE3 before
sum((outcomes2-beta0)^2)/length(outcomes2)
mean(outcomes5)- beta0   # MHDE3 after
sum((outcomes5-beta0 )^2)/length(outcomes5)

mean(outcomes1)- beta0   # MHDE4 before
sum((outcomes1-beta0)^2)/length(outcomes1)
mean(outcomes4)- beta0   # MHDE4 after
sum((outcomes4-beta0 )^2)/length(outcomes4)

mean(outcomes3)- beta0   # PMLE before
sum((outcomes3-beta0)^2)/length(outcomes3)
mean(outcomes6)- beta0   # PMLE after
sum((outcomes6-beta0)^2)/length(outcomes6)
proc.time()-ptm


