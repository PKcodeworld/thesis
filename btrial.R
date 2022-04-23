rm(list=ls())
library(KMsurv) 
library(survival)
data("btrial")


alldata = btrial
alldata$im=alldata$im-1
#rescale im = 1 if positive, im = 0 if negative


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

alldata0 = subset(alldata,im==0)   #group 0 (z=0)   
alldata1 = subset(alldata,im==1)   #group 1

library(bshazard)
fit0 = bshazard(Surv(time,death)~1,data=alldata0)
fit1 = bshazard(Surv(time,death)~1,data=alldata1)
baseline1 <-function(t) {kernelreg(fit0$time,fit0$hazard,t)}

Hellinger_distance1 <-function(beta){
  b = beta
  alltime = fit1$time
  all_0hazard = fit1$hazard/exp(1*b)
  
  sum = 0
  for (i in 1:length(alltime)){
    sum = sum + (baseline1(alltime[i])^(1/2)-all_0hazard[i]^(1/2))^2
  }
  sum
}


library(muhaz)
mfit0<-muhaz(alldata0$time,alldata0$death,min.time = min(alldata0$time),max.time = max(alldata0$time))
baseline0 <-function(t) {kernelreg(mfit0$est.grid,mfit0$haz.est,t)}
mfit1 = muhaz(alldata1$time,alldata1$death,min.time = min(alldata1$time),max.time = max(alldata1$time),n.est.grid=length(alldata1$time))
haz1<-function(t) {kernelreg(mfit1$est.grid,mfit1$haz.est,t)}

Hellinger_distance0 <-function(beta){
  b = beta
  alltime = fit1$time
  all_0hazard = haz1(alldata1$time)/exp(1*b)
  
  sum = 0
  for (i in 1:length(alltime)){
    sum = sum + (baseline0(alltime[i])^(1/2)-all_0hazard[i]^(1/2))^2
  }
  sum
}
hat0 = optimize(Hellinger_distance0,lower = 0,upper = 5)
hat0$minimum  # muhaz (non-parametric kernel estimate of conditional hazard)

hat1 = optimize(Hellinger_distance1,lower = 0,upper = 10)
hat1$minimum  #bashazard (parametric estimate of conditional hazard)

coxfit = coxph(Surv(time,death)~im,data=alldata)
coxfit$coefficients[1]
summary(coxfit)
