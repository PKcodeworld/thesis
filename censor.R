############# initial setting ########
beta0 = 0.8
#h0(t)=0.1+0.2t
#Z = 0,1,2
rm(list=ls())

#____________generate data_____________________
N=300   # 50 - 100 - 300
cen_rate = 0.3
# censoring rate 10%-30%
# repeat 500 times calculate MSE/Bias
inversion_gen<-function(n,Z){
  temp = runif(n)
  x = (-1+sqrt(1-40*log(1-temp)/exp(beta0*Z)))/2
}

library(bshazard)

Hellinger_distance <-function(beta){
  b = beta
  alltime = c(fit0$time,fit1$time,fit2$time)
  all_0hazard=c(fit0$hazard/exp(0*b),fit1$hazard/exp(1*b),fit2$hazard/exp(2*b))
  
  hazfun <- function(t){
    kernelreg(alltime,all_0hazard,t)
  }
  
  distance0 <- function(x){
    (sqrt(haz0fun(x))-sqrt(hazfun(x)))^2
  }
  
  distance1 <- function(x){
    (sqrt(haz1fun(x))-sqrt(hazfun(x)*exp(b*1)))^2
  }
  
  
  distance2 <- function(x){
    (sqrt(haz2fun(x))-sqrt(hazfun(x)*exp(b*2)))^2
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
  
  up0 = max(newdata0)
  up1 = max(newdata1)
  up2 = max(newdata2)
  integrate(distance0,low0,up0)$value+integrate(distance1,low1,up1)$value+integrate(distance2,low2,up2)$value
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

outcomes = c()

for (i in 1:500){
newdata0 = inversion_gen(N,0)  # 0-8
for (i in 1:(N*cen_rate)){
  temp = 100
  while (temp>newdata0[i]){
    temp = runif(1,0,10)
  }
  newdata0[i] = temp
}

newdata1 = inversion_gen(N,1)  # 0-6
for (i in 1:(N*cen_rate)){
  temp = 100
  while (temp>newdata1[i]){
    temp = runif(1,0,10)
  }
  newdata1[i] = temp
}
newdata2 = inversion_gen(N,2)  # 0-4
for (i in 1:(N*cen_rate)){
  temp = 100
  while (temp>newdata2[i]){
    temp = runif(1,0,10)
  }
  newdata2[i] = temp
}


del = c(rep(0,N*cen_rate),rep(1,N*(1-cen_rate)) )   #censoring indicator
fit0<-bshazard(Surv(newdata0,del)~1)
fit1<-bshazard(Surv(newdata1,del)~1)
fit2<-bshazard(Surv(newdata2,del)~1)

#########################
haz0fun <-function(t) {kernelreg(fit0$time,fit0$hazard,t)}
haz1fun <-function(t) {kernelreg(fit1$time,fit1$hazard,t)}
haz2fun <-function(t) {kernelreg(fit2$time,fit2$hazard,t)}

#b = seq(0,0.7,0.01)
#dis = rep(0,length(b))
#for (i in 1:length(b)){
#  dis[i] = Hellinger_distance(b[i])
#}
#plot(b,dis)

hat = optimize(Hellinger_distance,lower = 0,upper = 3)
outcomes = c(outcomes,hat$minimum)
}

mean(outcomes)- beta0
sum((outcomes-mean(outcomes))^2)/length(outcomes)
