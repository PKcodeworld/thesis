############# initial setting ########
beta1 = 0.2
beta2 = 0.8
#h0(t)=0.1+0.2t
#Z1 = 0,1, Z2=1,3
rm(list=ls())
#############  test area #############
p <- function(x){
  (0.1+0.2*x)*exp(-(0.1*x+0.1*x^2))
}
p2<-function(x){
  exp(0.5)*(0.1+0.2*x)*exp(-exp(0.5)*(0.1*x+0.1*x^2))
}

hist(newdata1)

x = seq(0,8,0.01)

par(mfrow=c(1,1))
plot(x,p2(x))
p2(0)

xx = seq(0,8,0.01)
plot(xx,p(xx))

plot(fit$time,fit$hazard)
max(fit$time)

rm(list=ls())
#____________generate data_____________________
N=10000   # 50 - 100 - 300
# censoring rate 10%-30%
# repeat 500 times calculate MSE/Bias
inversion_gen<-function(n,Z1,Z2){
  temp = runif(n)
  x = (-1+sqrt(1-40*log(1-temp)/exp(beta1*Z1+beta2*Z2)))/2
}

newdata0 = inversion_gen(N,0,1)  
newdata1 = inversion_gen(N,0,3)  
newdata2 = inversion_gen(N,1,1)  
newdata3 = inversion_gen(N,1,3)
#______hazard estimation based on interval-Poisson assumption_______

library(bshazard)

del = rep(1,N)    #censoring indicator
fit0<-bshazard(Surv(newdata0,del)~1)
fit1<-bshazard(Surv(newdata1,del)~1)
fit2<-bshazard(Surv(newdata2,del)~1)
fit3<-bshazard(Surv(newdata3,del)~1)

#########################
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
#########################
haz0fun <-function(t) {kernelreg(fit0$time,fit0$hazard,t)}
haz1fun <-function(t) {kernelreg(fit1$time,fit1$hazard,t)}
haz2fun <-function(t) {kernelreg(fit2$time,fit2$hazard,t)}
haz3fun <-function(t) {kernelreg(fit3$time,fit3$hazard,t)}

Hellinger_distance <-function(beta){
  b = beta
  
  alltime = c(fit0$time,fit1$time,fit2$time,fit3$time)
  all_0hazard=c(fit0$hazard/exp(b[2]),fit1$hazard/exp(3*b[2]),fit2$hazard/exp(b[1]+b[2]),fit3$hazard/exp(b[1]+3*b[2]))
  
  
  
  hazfun <- function(t){
    kernelreg(alltime,all_0hazard,t)
  }
  
  distance0 <- function(x){
    (sqrt(haz0fun(x)/exp(b[2]))-sqrt(hazfun(x)))^2
  }
  
  distance1 <- function(x){
    (sqrt(haz1fun(x)/exp(3*b[2]))-sqrt(hazfun(x)))^2
  }
  
  
  distance2 <- function(x){
    (sqrt(haz2fun(x)/exp(b[1]+b[2]))-sqrt(hazfun(x)))^2
  }
  
  distance3 <- function(x){
    (sqrt(haz2fun(x)/exp(b[1]+3*b[2]))-sqrt(hazfun(x)))^2
  }
  
  # x=seq(0,10,0.1)
  # y=rep(0,length(x))
  # for (i in 1:length(x)){
  #  y[i]=distance0(x[i])
  # }
  # plot(x,y)
  low0=min(newdata0)
  low1=min(newdata1)
  low2=min(newdata2)
  low3=min(newdata3)
  
  up0 = max(newdata0)
  up1 = max(newdata1)
  up2 = max(newdata2)
  up3 = max(newdata3)
  integrate(distance0,low0,up0)$value+integrate(distance1,low1,up1)$value+integrate(distance2,low2,up2)$value+integrate(distance3,low3,up3)$value
}



#########################################

alltime = c(fit0$time,fit1$time,fit2$time,fit3$time)
alldel = rep(1,4*N)
g1 = c(rep(0,N),rep(0,N),rep(1,N),rep(1,N))
g2 = c(rep(1,N),rep(3,N),rep(1,N),rep(3,N))
library(survival)
coxfit = coxph(Surv(alltime,alldel)~g1+g2)
coxfit$coefficients

nlm(Hellinger_distance,c(0.2,0.2))
