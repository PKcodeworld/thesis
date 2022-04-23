rm(list=ls())

beta0 = 1
N=100  
library(survival)
hx=N^(-1/5)
ht=N^(-1/5)
scale_size1 = 2
scale_size2 = 5

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



Pgood = c()
Hgood = c()
Pscale = c()
Hscale1 = c()
Hscale2 = c()

ptm<-proc.time()
for (i in 1:100)
  
{
  Z = runif(N,0.1,2)
  alltime= c()
  for (i in 1:N){
    temp = inversion_gen(1,Z[i])
    alltime = c(alltime,temp)
  }
  
  del = rep(1,N)
   hat1=optimize(MHD,lower = 0,upper = 10)
   H_origin_beta = hat1$minimum
   Hgood = c(Hgood,H_origin_beta)
  
  coxfit = coxph(Surv(alltime,del)~Z)
  P_origin_beta=coxfit$coefficients[1]
  Pgood = c(Pgood,P_origin_beta)
  
  Z = Z*scale_size1
  
  coxfit = coxph(Surv(alltime,del)~Z)
  P_scale_beta=coxfit$coefficients[1]
  
  Pscale = c(Pscale,P_scale_beta)  
  
   hat2=optimize(MHD,lower = 0,upper = 10)
   H_scale_beta=hat2$minimum
   Hscale1 = c(Hscale1,H_scale_beta)
   
   Z=Z/scale_size1*scale_size2
   
   hat3=optimize(MHD,lower = 0,upper = 10)
   H_scale_beta=hat3$minimum
   Hscale2 = c(Hscale2,H_scale_beta)
}


proc.time()-ptm

mean(Pgood)
var(Pgood)
mean(Pscale)*scale_size1
var(Pscale)

mean(Hgood)
var(Hgood)
mean(Hscale1)*scale_size1
var(Hscale1)

mean(Hscale2)*scale_size2
var(Hscale2)
