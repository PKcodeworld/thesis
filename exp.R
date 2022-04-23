rm(list=ls())
############# initial setting ########
beta0 = 0.5
beta1 = 0.5
#h0(t)=0.1+0.2t
#Z0 = runif(N,0,2)
#Z1 = runif(N,0,2)


#____________generate data_____________________
N=300   # 50 - 100 - 300
# censoring rate 10%-30%
# repeat 500 times calculate MSE/Bias
inversion_gen<-function(n,Z0,Z1){
  temp = runif(n)
  x = (-1+sqrt(1-40*log(1-temp)/exp(beta0*Z0+beta1*Z1)))/2
}



condional_haz = function(Covar1,Covar2,times,del,t,x1,x2){
  
  
  hx=N^(-1/5)
  ht=N^(-1/5)
  
  n= length(times)
  K = function(x) {1/sqrt(2*pi) * exp(-sum(x^2)/2)}
  K1 = function(x){1/sqrt(2*pi) * exp(-x^2/2)}
  W = function(y) { integrate(K1,lower=-Inf,upper=y)$value }
  
  #The optimal bandwidths are both of order n^(???1/6)
  #and the corresponding AMISE is of order n^(???2/3)
  
  
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


#______hazard estimation based on interval-Poisson assumption_______

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

baseline <- function(t){0.1+0.2*t}
#baseline = haz0fun

Hellinger_distance <-function(beta){
  
  b0 = beta[1]
  b1 = beta[2]
  
  
  expset = rep(0,N)
  estset = rep(0,N)
  for (i in 1:N){
    haz= condional_haz(Z0,Z1,alltime,alldel,alltime[i],Z0[i],Z1[i])
    expset[i]=haz/ baseline(alltime[i])
    estset[i]=exp(b0*Z0[i]+b1*Z1[i])
  }
  
  # est_haz = function(t){kernelreg(alltime,zeroset,t)}
  
  
  sum = 0
  for (i in 1:length(alltime)){
    sum = sum + (expset[i]^(1/2)-estset[i]^(1/2))^2
  }
  sum
}

########################################
#PMLE works very well

b0 = c()
b1 = c()
Hb0 = c()
Hb1 = c()

ptm=proc.time()
for(i in 1:100){
  
  Z0 = runif(N,0,1)
  Z1 = runif(N,0,1)
  
  alltime = rep(0,N)
  for (i in 1:N){
    alltime[i] = inversion_gen(1,Z0[i],Z1[i]) 
  }
  
  del = rep(1,N)    #censoring indicator
  alldel = rep(1,N)
  es=optim(c(0,0),Hellinger_distance)
  es$par
  Hb0 = c(Hb0,es$par[1])
  Hb1 = c(Hb1,es$par[2])
  
  coxfit = coxph(Surv(alltime,alldel)~Z0+Z1)
  b0 = c(b0,coxfit$coefficients[1])
  b1 = c(b1,coxfit$coefficients[2])
  
}

proc.time()-ptm
mean(Hb0)- beta0
sum((Hb0-mean(Hb0))^2)/length(Hb0)
mean(Hb1) - beta1
sum((Hb1-mean(Hb1))^2)/length(Hb1)

mean(b0)- beta0
sum((b0-mean(b0))^2)/length(b0)
mean(b1) - beta1
sum((b1-mean(b1))^2)/length(b1)

