rm(list=ls())
library(KMsurv)
library(survival)
data(hodg)
burn = hodg

#gtype: Graft type (1=allogenic, 2=autologous) 
#dtype: Disease type (1=Non Hodgkin lymphoma, 2=Hodgkins disease) 
#time: Time to death or relapse, days 
#delta: Death/relapse indicator (0=alive, 1=dead) 

#rescale the covariates
burn$time = burn$time
burn$gtype = burn$gtype-1  #gtype: Graft type (0=allogenic, 1=autologous) 
burn$dtype = burn$dtype-1 #dtype: Disease type (0=Non Hodgkin lymphoma, 1=Hodgkins disease) 
data0 = subset(burn,gtype ==0&dtype==0)
data1 = subset(burn,gtype ==0&dtype==1)
data2 = subset(burn,gtype ==1&dtype==0)
data3 = subset(burn,gtype ==1&dtype==1)


#kernel regression  ##############################################
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

######################################################

library(muhaz)
fit0<-muhaz(data0$time,data0$delta,min.time=min(burn$time),max.time=max(burn$time))
plot(fit0$est.grid,fit0$haz.est)
baseline <-function(t) {kernelreg(fit0$est.grid,fit0$haz.est,t)}

xx=seq(0,2000,1)
plot(xx,baseline(xx))
fit1<-muhaz(data1$time,data1$delta,min.time = min(data1$time),max.time = max(data1$time))
haz1<-function(t) {kernelreg(fit1$est.grid,fit1$haz.est,t)}
fit2<-muhaz(data2$time,data2$delta,min.time = min(data2$time),max.time = max(data2$time))
haz2<-function(t) {kernelreg(fit2$est.grid,fit2$haz.est,t)}
fit3<-muhaz(data3$time,data3$delta,min.time = min(data3$time),max.time = max(data3$time))
haz3<-function(t) {kernelreg(fit3$est.grid,fit3$haz.est,t)}
#hazard0 = baseline(data0$time)
hazard1 = haz1(data1$time)
plot(data3$time,hazard3)
hazard2 = haz2(data2$time)
hazard3 = haz3(data3$time)

Hellinger_distance <-function(beta){

  b0 = beta[1]
  b1 = beta[2]
  
 # alltime = burn$time
  alltime = c(data1$time,data2$time,data3$time)
 # all_0hazard=c(hazard0,hazard1/exp(b1),hazard2/exp(b0),hazard3/exp(b0+b1))
  all_0hazard=c(hazard1/exp(b1),hazard2/exp(b0),hazard3/exp(b0+b1))
  sum = 0
  for (i in 1:length(alltime)){
    sum = sum + (baseline(alltime[i])^(1/2)-all_0hazard[i]^(1/2))^2
  }
  sum
}
es=optim(c(-1,1),Hellinger_distance)
es$par
coxfit = coxph(Surv(time,delta)~gtype+dtype,data=burn)
summary(coxfit)
