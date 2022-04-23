############# initial setting ########
beta0 = 0.5
beta1 = 0.8
#h0(t)=0.1+0.2t
#Z0 = 0,1
#Z1 = 0,1
rm(list=ls())

#____________generate data_____________________
N=50   # 50 - 100 - 300
# censoring rate 10%-30%
# repeat 500 times calculate MSE/Bias
inversion_gen<-function(n,Z0,Z1){
  temp = runif(n)
  x = (-1+sqrt(1-40*log(1-temp)/exp(beta0*Z0+beta1*Z1)))/2
}

newdata0 = inversion_gen(N,0,0)  
newdata01 = inversion_gen(N,0,1)  
newdata10 = inversion_gen(N,1,0)
newdata11 = inversion_gen(N,1,1)

#______hazard estimation based on interval-Poisson assumption_______

library(bshazard)

del = rep(1,N)    #censoring indicator
fit0<-bshazard(Surv(newdata0,del)~1)
fit01<-bshazard(Surv(newdata01,del)~1)
fit10<-bshazard(Surv(newdata10,del)~1)
fit11<-bshazard(Surv(newdata11,del)~1)

par(mfrow=c(1,1))
plot(fit0$time,fit0$hazard)
plot(fit01$time,fit01$hazard)
plot(fit10$time,fit10$hazard)
plot(fit11$time,fit11$hazard)

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
#baseline <- function(t){0.1+0.2*t}
baseline = haz0fun

Hellinger_distance <-function(beta){
  b0 = beta[1]
  b1 = beta[2]
  alltime = c(fit0$time,fit01$time,fit10$time,fit11$time)
  all_0hazard=c(fit0$hazard,fit01$hazard/exp(b1),fit10$hazard/exp(b0),fit11$hazard/exp(b0+b1))
  
  sum = 0
  for (i in 1:length(alltime)){
    sum = sum + (baseline(alltime[i])^(1/2)-all_0hazard[i]^(1/2))^2
  }
  sum
}


Hb0 = c()
Hb1 = c()
for (i in 1:300){
  newdata0 = inversion_gen(N,0,0)  
  newdata01 = inversion_gen(N,0,1)  
  newdata10 = inversion_gen(N,1,0)
  newdata11 = inversion_gen(N,1,1)
  
  del = rep(1,N)    #censoring indicator
  fit0<-bshazard(Surv(newdata0,del)~1)
  fit01<-bshazard(Surv(newdata01,del)~1)
  fit10<-bshazard(Surv(newdata10,del)~1)
  fit11<-bshazard(Surv(newdata11,del)~1)

es = nlm(Hellinger_distance,c(0.6,0.6))

Hb0 = c(Hb0,es$estimate[1])
Hb1 = c(Hb1,es$estimate[2])

}
mean(Hb0)-beta0
mean(Hb1)-beta1
sum((Hb0-mean(Hb0))^2)/length(Hb0)
sum((Hb1-mean(Hb1))^2)/length(Hb1)


########################################
#PMLE works very well

b0 = c()
b1 = c()

for(i in 1:300){
  
newdata0 = inversion_gen(N,0,0)  
newdata01 = inversion_gen(N,0,1)  
newdata10 = inversion_gen(N,1,0)
newdata11 = inversion_gen(N,1,1)

del = rep(1,N)    #censoring indicator
fit0<-bshazard(Surv(newdata0,del)~1)
fit01<-bshazard(Surv(newdata01,del)~1)
fit10<-bshazard(Surv(newdata10,del)~1)
fit11<-bshazard(Surv(newdata11,del)~1)

es = nlm(Hellinger_distance,c(0.6,0.6))

Hb0 = c(Hb0,es$estimate[1])
Hb1 = c(Hb1,es$estimate[2])

alltime = c(newdata0,newdata01,newdata10,newdata11)
alldel = rep(1,4*N)
g0 = c(rep(0,N),rep(0,N),rep(1,N),rep(1,N))
g1 = c(rep(0,N),rep(1,N),rep(0,N),rep(1,N))
library(survival)
coxfit = coxph(Surv(alltime,alldel)~g0+g1)
b0 = c(b0,coxfit$coefficients[1])
b1 = c(b1,coxfit$coefficients[2])
}

mean(b0)- beta0
mean(b1) - beta1
sum((b0-mean(b0))^2)/length(b0)
sum((b1-mean(b1))^2)/length(b1)

