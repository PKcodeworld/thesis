rm(list=ls())

beta0 = 0.5
N=100 
library(survival)

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

outcomes1 = c()
outcomes2 = c()
outcomes3 = c()
outcomes4 = c()
ptm<-proc.time()
for (i in 1:100){
  
  Z = runif(N,0,4)
  alltime= c()
  for (i in 1:N){
    temp = inversion_gen(1,Z[i])
    alltime = c(alltime,temp)
  }
  summary(alltime)
  del=rep(1,N)
  # 20% censor
  for (i in 1:N){
    temp = runif(1,0,4.2)
    if (temp<alltime[i]){
      alltime[i]=temp
      del[i]=0
    }
  }

  
  hx=N^(-1/5)
  ht=N^(-1/5)
  
  hat=optimize(MHD1,lower = 0,upper = 10)
  outcomes1 = c(outcomes1 , hat$minimum)
  
  
  hat=optimize(MHD3,lower = 0,upper = 10)
  outcomes3 = c(outcomes3 , hat$minimum)
  
  
  fit = coxph(Surv(alltime,del)~Z)
  outcomes4 = c(outcomes4,fit$coefficients[1])
  
} 
mean(outcomes1)- beta0   # cons_h0 vs est_h0  
sum((outcomes1-beta0)^2)/length(outcomes1)

#mean(outcomes2)- beta0   # zeroset vs est_h0
#sum((outcomes2-beta0)^2)/length(outcomes2)

mean(outcomes3)- beta0   # zeroset vs  known
sum((outcomes3-beta0 )^2)/length(outcomes3)

mean(outcomes4)- beta0   # PMLE 
sum((outcomes4-beta0)^2)/length(outcomes4)
proc.time()-ptm
