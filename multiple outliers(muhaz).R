############# initial setting ########

#h0(t)=0.1+0.2t
#Z = 0,1,2
rm(list=ls())
beta0 = 1

#____________generate data_____________________
N=300   # 50 - 100 - 300

inversion_gen<-function(n,Z){
  temp = runif(n)
  x = (-1+sqrt(1-40*log(1-temp)/exp(beta0*Z)))/2
}


Ybar<-function(x,dataset){
  sum(1*(dataset>=x))/length(dataset)
}

#______hazard estimation based on interval-Poisson assumption_______

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

#baseline <-function(t) 0.1+0.2*t


Hellinger_distance0 <-function(beta){
  b = beta
  alltime = c(fit0$est.grid,fit1$est.grid,fit2$est.grid)
  all_0hazard=c(fit0$haz.est/exp(0*b),fit1$haz.est/exp(1*b),fit2$haz.est/exp(2*b))
  baseline  <-function(t) {kernelreg(alltime,all_0hazard,t)}
  sum = 0
  for (i in 1:length(alltime)){
    sum = sum + (baseline(alltime[i])^(1/2)-all_0hazard[i]^(1/2))^2
  }
  sum
}
max(all_0hazard)
xx= seq(0,50,0.1)
plot(xx,baseline(xx))

outcomes1 = c()
outcomes2 = c()
outcomes3 = c()
outcomes4 = c()

for (i in 1:100){
  
  newdata0 = inversion_gen(N,0)  # 0-8
  newdata1 = inversion_gen(N,1)  # 0-6
  newdata2 = inversion_gen(N,2)  # 0-4
  
  del = rep(1,N)    #censoring indicator
  fit0<-muhaz(newdata0,del,n.est.grid = N)
  fit1<-muhaz(newdata1,del,n.est.grid = N)
  fit2<-muhaz(newdata2,del,n.est.grid = N)
  
  
  hat = optimize(Hellinger_distance0,lower = 0,upper = 5)
  outcomes1 = c(outcomes1,hat$minimum)
  
  alltime = c(newdata0,newdata1,newdata2)
  alldel = rep(1,(3*N))
  g = c(rep(0,N),rep(1,N),rep(2,N))
  
  coxfit = coxph(Surv(alltime,alldel)~g)
  outcomes2 = c(outcomes2,coxfit$coefficients[1])
  
  
  # setting outliners
  
  percent=0.05
  # setting outliners
  for (i in 1:floor(N*percent)){
    newdata2[i]= rnorm(1,40,1)
  }
  
  fit2<-muhaz(newdata2,del,n.est.grid = N)
  
  hat = optimize(Hellinger_distance0,lower = 0,upper = 5)
  outcomes3 = c(outcomes3,hat$minimum)
  
  alltime = c(newdata0,newdata1,newdata2)
  coxfit = coxph(Surv(alltime,alldel)~g)
  
  outcomes4 = c(outcomes4,coxfit$coefficients[1])
  
}

mean(outcomes1)- beta0   #new before
sum((outcomes1-mean(outcomes1))^2)/length(outcomes1)

mean(outcomes3)- beta0    #new after
sum((outcomes3-mean(outcomes3))^2)/length(outcomes3)

mean(outcomes2)- beta0    #PMLE before
sum((outcomes2-mean(outcomes2))^2)/length(outcomes2)

mean(outcomes4)- beta0    #PMLE after
sum((outcomes4-mean(outcomes4))^2)/length(outcomes4)
