rm(list=ls())

beta0 = 0.5
N=50  
library(survival)
hx=N^(-1/5)
ht=N^(-1/5)

#par(mfrow=c(1,2))

condional_haz = function(Covar,times,del,t,x){
  
  n= length(times)
  K = function(x) {1/sqrt(2*pi) * exp(-(x^2)/2)}
  W = function(y) { integrate(K,lower=-Inf,upper=y)$value  }
  
  #The optimal bandwidths are both of order n^(???1/6)
  #and the corresponding AMISE is of order n^(???2/3)
  
  
  wx=rep(0,n)
  sumw = 0
  for (i in 1:n){
    sumw = sumw+K((x-Covar[i])/hx)
  }
  sumdeno = 0
  sumnumer = 0
  
  for (i in 1:n){
    wx[i] = K((x-Covar[i])/hx)/sumw
    sumdeno = sumdeno+wx[i]*W((times[i]-t)/ht)
    sumnumer = sumnumer+wx[i]*K((t-times[i])/ht)*del[i]
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



inversion_gen<-function(n,Z){
  temp = runif(n)
  x = (-1+sqrt(1-40*log(1-temp)/exp(beta0*Z)))/2
}

MHD = function(b){
  
  zeroset= rep(0,length(alltime))
  for (i in 1:length(alltime)){
    zeroset[i]=condional_haz(Z,alltime,del,alltime[i],Z[i])/exp(b*Z[i])
  }
  
  cons_h0 = function(t){kernelreg(alltime,zeroset,t)}
  
  #est_h0= function(t){condional_haz(Z,alltime,del,t,0)}
  est_h0 = function(t) 0.1+0.2*t
  distance= 0
  for (i in 1:length(alltime)){
    distance = distance+(cons_h0(alltime[i])^0.5-est_h0(alltime[i])^0.5)^2
  }
  distance
}


P0 = c()
P20 = c()
P40 = c()
H0 = c()
H20 = c()
H40 = c()

P0out = c()
P20out = c()
P40out = c()
H0out = c()
H20out = c()
H40out = c()

rate20 = c()
rate40 = c()

ptm<-proc.time()
rate = c()
for (i in 1:100){
  Z = runif(N,0.2,2)
  
  alltime= c()
  for (i in 1:N){
    temp = inversion_gen(1,Z[i])
    alltime = c(alltime,temp)
  }
  del = rep(1,N)
  
  # backup
  alltime_ori = alltime
  del_ori=del
  Z_ori = Z  
  
  coxfit0 = coxph(Surv(alltime,del)~Z)
  P0 = c(P0,coxfit0$coefficients[1])
  
  hat0=optimize(MHD,lower = 0,upper = 10)
  H0 = c(H0,hat0$minimum)
  
  # one fixed outlier
  index = sample(which(del==1),1)
  Z[index]=5
  coxfit0out = coxph(Surv(alltime,del)~Z)
  P0out = c(P0out,coxfit0out$coefficients[1])
  
  hat0out=optimize(MHD,lower = 0,upper = 10)
  H0out = c(H0out,hat0out$minimum)
  
  #recover
  alltime = alltime_ori
  del=del_ori
  Z = Z_ori 
  #########################################
  #censor = 20%
  for (i in 1:N){
    temp = runif(1,0,8.4)
    if (temp<alltime[i]){
      alltime[i]=temp
      del[i]=0
    }
  }
  
  rate20 = c(rate20,1-sum(del)/N)
  coxfit20 = coxph(Surv(alltime,del)~Z)
  P20= c(P20,coxfit20$coefficients[1])
  
  hat20=optimize(MHD,lower = 0,upper = 10)
  H20 = c(H20,hat20$minimum)
  
  # one fixed outlier
  index = sample(which(del==1),1)
  Z[index]=5
  coxfit20out = coxph(Surv(alltime,del)~Z)
  P20out = c(P20out,coxfit20out$coefficients[1])
  
  hat20out=optimize(MHD,lower = 0,upper = 10)
  H20out = c(H20out,hat20out$minimum)
  
  #recover
  alltime = alltime_ori
  del=del_ori
  Z = Z_ori 
  #########################################
  #censor = 40%
  for (i in 1:N){
    temp = runif(1,0,4.2)
    if (temp<alltime[i]){
      alltime[i]=temp
      del[i]=0
    }
  }
  
  rate40 = c(rate40,1-sum(del)/N)
  coxfit40 = coxph(Surv(alltime,del)~Z)
  P40= c(P40,coxfit40$coefficients[1])
  
  hat40=optimize(MHD,lower = 0,upper = 10)
  H40 = c(H40,hat40$minimum)
  
  # one fixed outlier
  index = sample(which(del==1),1)
  Z[index]=5
  coxfit40out = coxph(Surv(alltime,del)~Z)
  P40out = c(P40out,coxfit40out$coefficients[1])
  
  hat40out=optimize(MHD,lower = 0,upper = 10)
  H40out = c(H40out,hat40out$minimum)
  
 
  
  ######################################
  
  
}
mean(rate20)
mean(rate40)
proc.time()-ptm

mean(P0)-beta0
sum((P0-mean(P0))^2)/N
mean(P20)-beta0
sum((P20-mean(P20))^2)/N
mean(P40)-beta0
sum((P40-mean(P40))^2)/N

mean(H0)-beta0
sum((H0-mean(H0))^2)/N
mean(H20)-beta0
sum((H20-mean(H20))^2)/N
mean(H40)-beta0
sum((H40-mean(H40))^2)/N

######################
mean(P0out)-beta0
sum((P0out-mean(P0out))^2)/N
mean(P20out)-beta0
sum((P20out-mean(P20out))^2)/N
mean(P40out)-beta0
sum((P40out-mean(P40out))^2)/N

mean(H0out)-beta0
sum((H0out-mean(H0out))^2)/N
mean(H20out)-beta0
sum((H20out-mean(H20out))^2)/N
mean(H40out)-beta0
sum((H40out-mean(H40out))^2)/N
