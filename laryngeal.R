rm(list=ls())
library(KMsurv)
library(survival)
data("larynx")
lary = larynx

#stage: Disease stage (1-4) from TNM cancer staging classification.
#time: Time from first treatment until death, or end of study.
#age: Age at diagnosis.
#delta: Indicator of death [1, if patient died at time t; 0, otherwise]. 

#rescale the covariates
lary$stage=lary$stage-1 #stage: Disease stage (0-3)
lary$age = lary$age-min(lary$age) #start from 0
coxfit = coxph(Surv(time,delta)~stage+age,data=lary)
summary(coxfit)

#############################################################

condional_haz = function(Covar1,Covar2,times,del,t,x1,x2){
  
  hx=N^(-1/5)
  ht=N^(-1/5)
  
  n= length(times)
  K = function(x) {1/sqrt(2*pi) * exp(-sum(x^2)/2)}
  K1 = function(x){1/sqrt(2*pi) * exp(-x^2/2)}
  W = function(y) { integrate(K1,lower=-Inf,upper=y)$value }
  
  #The optimal bandwidths are both of order n^(1/6)
  #and the corresponding AMISE is of order n^(2/3)
  
  wx=rep(1,n)
  sumw = 0
  for (i in 1:n){
  
    sumw = sumw+K(c((x1-Covar1[i])/hx,(x2-Covar2[i])/hx))
  }

  sumdeno = 0
  sumnumer = 0
  for (i in 1:n){
    
    wx[i] = K(c((x1-Covar1[i])/hx,(x2-Covar2[i])/hx))/sumw
    sumdeno = sumdeno+wx[i]*W((times[i]-t)/ht)
    sumnumer = sumnumer+wx[i]*K1((t-times[i])/ht)*del[i]
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
  
N = nrow(lary)
hazard = rep(0,N)
zeroset = rep(0,N)
for (i in 1:N){
    hazard[i] = condional_haz(lary$stage,lary$age,lary$time,lary$delta,lary$time[i],lary$stage[i],lary$age[i])
}

baseline = function(t){condional_haz(lary$stage,lary$age,lary$time,lary$delta,t,0,0)}

Hellinger_distance <-function(beta){
  b0 = beta[1]
  b1 = beta[2]

  for (i in 1:N){
    zeroset[i]=hazard[i]/exp(b0*lary$stage[i]+b1*lary$age[i])
  }
  
  sum = 0
  for (i in 1:N){
    sum = sum + (baseline(lary$time[i])^(1/2)-zeroset[i]^(1/2))^2
  }
  sum
}
es=optim(c(0,0),Hellinger_distance)
es$par

