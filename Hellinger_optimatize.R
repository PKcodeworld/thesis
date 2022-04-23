rm(list=ls())
beta0 = 2

#____________generate data_____________________
N=50   # 50 - 100 - 300

inversion_gen<-function(n,Z){
  temp = runif(n)
  x = (-1+sqrt(1-40*log(1-temp)/exp(beta0*Z)))/2
}

#newdata0 = inversion_gen(N,0)  # 0-8
#newdata1 = inversion_gen(N,1)  # 0-6
#newdata2 = inversion_gen(N,2)  # 0-4

#del = rep(1,N)    #censoring indicator


library(bshazard)
library(survival)
library(muhaz)


#########################
kernelreg<-function(x,y,t){
  n <- length(x)
  K <- function(x)  1/sqrt(2*pi) * exp(-(x^2)/2)
  h= n^(-0.2)  #Silverman's Rule of Thumb
  numer=0
  deno = 0
  for(i in 1:n){
    numer = numer+y[i]*K((x[i]-t)/h)
    deno= deno+K((x[i]-t)/h)
  }
  return (numer/deno)
  
}



#########################
haz0fun <-function(t) {kernelreg(fit0$est.grid,fit0$haz.est,t)}
haz1fun <-function(t) {kernelreg(fit1$est.grid,fit1$haz.est,t)}
haz2fun <-function(t) {kernelreg(fit2$est.grid,fit2$haz.est,t)}



Hellinger_optimize <-function(init){
  b = init
  Distance <-function(beta){
    haz0fun <-function(t) {kernelreg(fit0$est.grid,fit0$haz.est,t)}
    haz1fun <-function(t) {kernelreg(fit1$est.grid,fit1$haz.est,t)}
    haz2fun <-function(t) {kernelreg(fit2$est.grid,fit2$haz.est,t)}
    baseline  <-function(t) {kernelreg(alltime,all_0hazard,t)}
    
    sum = 0
    for (i in 1:length(newdata1)){
      sum = sum + (baseline(newdata1[i])^(1/2)-(haz1fun(newdata1[i])/exp(1*beta))^(1/2))^2
    }
    for (i in 1:length(newdata2)){
      sum = sum + (baseline(newdata2[i])^(1/2)-(haz2fun(newdata2[i])/exp(2*beta))^(1/2))^2
    }
    sum
  }
  
 # estimation = c()
  
  for (i in 1:10){
  all_0hazard=c(fit0$haz.est/exp(0*b),fit1$haz.est/exp(1*b),fit2$haz.est/exp(2*b))
  
  h = optimize(Distance,lower = 0,upper = 10)
  b = h$minimum
  }
  b
}

outcomes1 = c()
outcomes2 = c()
ptm<-proc.time()
for (i in 1:100){
  
  newdata0 = inversion_gen(N,0)  # 0-8
  newdata1 = inversion_gen(N,1)  # 0-6
  newdata2 = inversion_gen(N,2)  # 0-4
  
  alltime = c(newdata0,newdata1,newdata2)
  alldel = rep(1,(3*N))
  g = c(rep(0,N),rep(1,N),rep(2,N))
  
  coxfit = coxph(Surv(alltime,alldel)~g)
  outcomes1 = c(outcomes1,coxfit$coefficients[1])
  
  del = rep(1,N)    #censoring indicator
  fit0<-muhaz(newdata0,del,n.est.grid = N)
  fit1<-muhaz(newdata1,del,n.est.grid = N)
  fit2<-muhaz(newdata2,del,n.est.grid = N)
  alltime = c(fit0$est.grid,fit1$est.grid,fit2$est.grid)


  Hb = Hellinger_optimize(5)
    outcomes2 = c(outcomes2,Hb)
  

  
  
  # setting outliners
  
 # percent=1/N
  # setting outliners
  #for (i in 1:floor(N*percent)){
 #   newdata2[i]= rnorm(1,40,1)
 # }
  
 # fit2<-muhaz(newdata2,del,n.est.grid = N)
  
 # hat = optimize(Hellinger_distance0,lower = 0,upper = 5)
 # outcomes3 = c(outcomes3,hat$minimum)
  
 # alltime = c(newdata0,newdata1,newdata2)
 # coxfit = coxph(Surv(alltime,alldel)~g)
  
 # outcomes4 = c(outcomes4,coxfit$coefficients[1])
  
}

mean(outcomes1)- beta0   #PMLE
sum((outcomes1-mean(outcomes1))^2)/length(outcomes1)


mean(outcomes2)- beta0    #new iteration
sum((outcomes2-mean(outcomes2))^2)/length(outcomes2)
proc.time() - ptm