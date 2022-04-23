rm(list=ls())

beta0 = 0.5
N=300  
library(survival)

condional_haz = function(Covar,times,del,t,x){
  
  n= length(times)
  K = function(x) {1/sqrt(2*pi) * exp(-(x^2)/2)}
  W = function(y) { integrate(K,lower=-Inf,upper=y)$value  }
  
  #The optimal bandwidths are both of order n^(-1/5)
  
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


MHD = function(b){
  zeroset= rep(0,length(alltime))
  est = rep(0,length(alltime))
  for (i in 1:length(alltime)){
    zeroset[i]=condional_haz(Z,alltime,del,alltime[i],Z[i])/exp(b*Z[i])
    est[i] = condional_haz(Z,alltime,del,alltime[i],0)
  }

  distance= 0
  for (i in 1:length(alltime)){
    distance = distance+(zeroset[i]^0.5-est[i]^0.5)^2
  }
  distance
}


outcomes1 = c()
outcomes2 = c()

ptm<-proc.time()
for (i in 1:100){
  Z = runif(N,0.1,2)
  alltime= c()
  for (i in 1:N){
    temp = inversion_gen(1,Z[i])
    alltime = c(alltime,temp)
  }
  del = rep(1,N)
  
  hx=N^(-1/5)*sd(Z)
  ht=N^(-1/5)*sd(alltime)
  hat=optimize(MHD,lower = 0,upper = 10)
  outcomes1 = c(outcomes1 , hat$minimum)
  
  fit = coxph(Surv(alltime,del)~Z)
  outcomes2 = c(outcomes2,fit$coefficients[1])
  
} 
mean(outcomes1)- beta0   # MHD
sum((outcomes1-beta0)^2)/length(outcomes1)

mean(outcomes2)- beta0   # PMLE 
sum((outcomes2-beta0)^2)/length(outcomes2)
proc.time()-ptm

