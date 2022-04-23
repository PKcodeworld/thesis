rm(list=ls())

############# initial setting ########
beta0 = 1
#h0(t)=0.1+0.2t
#Z = 0,1,2



#____________generate data_____________________
N=50   # 50 - 100 - 300



Ybar<-function(x,dataset){
  sum(1*(dataset>=x))/length(dataset)
}

inversion_gen<-function(n,Z){
  temp = runif(n)
  x = (-1+sqrt(1-40*log(1-temp)/exp(beta0*Z)))/2
}

library(bshazard)
library(survival)




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
baseline<-haz0fun

Hellinger_distance0 <-function(beta){
  b = beta
  alltime = c(fit0$time,fit1$time,fit2$time)
  
  distance0 <- function(x){
    (sqrt(haz0fun(x))-sqrt(baseline(x)))^2*Ybar(x,newdata0)
  }
  
  distance1 <- function(x){
    (sqrt(haz1fun(x))-sqrt(baseline(x)*exp(b*1)))^2*Ybar(x,newdata1)
  }
  
  
  distance2 <- function(x){
    (sqrt(haz2fun(x))-sqrt(baseline(x)*exp(b*2)))^2*Ybar(x,newdata2)
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

Hellinger_distance <-function(beta){
  b = beta
  alltime = c(fit0$time,fit1$time,fit2$time)
  all_0hazard=c(fit0$hazard/exp(0*b),fit1$hazard/exp(1*b),fit2$hazard/exp(2*b))
  
  sum = 0
  for (i in 1:length(alltime)){
    sum = sum + (baseline(alltime[i])^(1/2)-all_0hazard[i]^(1/2))^2
  }
  sum
}


outcomes1 = c()
outcomes2 = c()
outcomes3 = c()

for (i in 1:100){
  newdata0 = inversion_gen((N*9/5),0)  # 0-8
  newdata1 = inversion_gen((N*3/5),1)  # 0-6
  newdata2 = inversion_gen((N*3/5),2)  # 0-4
  
  del0 = rep(1,(N*9/5))    #censoring indicator
  del1 = rep(1,(N*3/5))
  del2 = rep(1,(N*3/5))
  fit0<-bshazard(Surv(newdata0,del0)~1)
  fit1<-bshazard(Surv(newdata1,del1)~1)
  fit2<-bshazard(Surv(newdata2,del2)~1)
  
  hat = optimize(Hellinger_distance0,lower = 0,upper = 1)
  outcomes1 = c(outcomes1,hat$minimum)
  
  hat = optimize(Hellinger_distance,lower = 0,upper = 1)
  outcomes2 = c(outcomes2,hat$minimum)
  
  alltime = c(newdata0,newdata1,newdata2)
  alldel = rep(1,3*N)
  g = c(rep(0,(N*9/5)),rep(1,(N*3/5)),rep(2,(N*3/5)))
  
  coxfit = coxph(Surv(alltime,alldel)~g)
  outcomes3 = c(outcomes3,coxfit$coefficients[1])
}

mean(outcomes1)-beta0   #Ying
sum((outcomes1-mean(outcomes1))^2)/length(outcomes1)

mean(outcomes2)-beta0    # discrete
sum((outcomes2-mean(outcomes2))^2)/length(outcomes2)

mean(outcomes3)-beta0
sum((outcomes3-mean(outcomes3))^2)/length(outcomes3)
