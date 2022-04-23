rm(list=ls())

beta0 = 0.5
N=200  
library(survival)


par(mfrow=c(1,2))

condional_haz = function(Covar,times,del,t,x){
  
  n= length(times)
  K = function(x) {1/sqrt(2*pi) * exp(-(x^2)/2)}
  W = function(y) { integrate(K,lower=-Inf,upper=y)$value  }

  #The optimal bandwidths are both of order n^(−1/6)
  #and the corresponding AMISE is of order n^(−2/3)

  
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

#newdata0 = inversion_gen(N,0)  # 0-8
#newdata1 = inversion_gen(N,1)  # 0-6
#newdata2 = inversion_gen(N,2)  # 0-4
#alltime = c(newdata0,newdata1,newdata2)

#Z= c(rep(0,N),rep(1,N),rep(2,N))
#del = rep(1,3*N)

#{Z = runif(N,0.5,2)
#alltime= c()
#for (i in 1:N){
#  temp = inversion_gen(1,Z[i])
#  alltime = c(alltime,temp)
#}
#del = rep(1,N)

#con0haz= function(t){condional_haz(Z,alltime,del,t,0)}


#xx= seq(0,max(alltime),0.1)
#yy=rep(0,length(xx))
#for (i in 1:length(xx)){
 # yy[i]=con0haz(xx[i])
#}#
#plot(xx,yy)



#zeroset= rep(0,length(alltime))
#for (i in 1:length(alltime)){
 # zeroset[i]=condional_haz(Z,alltime,del,alltime[i],Z[i])/exp(beta0*Z[i])
#}

#plot(alltime,zeroset)

#cons_h0 = function(t){kernelreg(alltime,zeroset,t)}
#plot(xx,cons_h0(xx))}

outcomes1 = c()
outcomes2 = c()
outcomes3 = c()
ptm<-proc.time()
for (i in 1:100){
Z = runif(N,0.5,4)
alltime= c()
for (i in 1:N){
  temp = inversion_gen(1,Z[i])
  alltime = c(alltime,temp)
}
del = rep(1,N)

hx=N^(-1/6)
ht=N^(-1/6)
hat=optimize(MHD,lower = 0,upper = 10)
outcomes1 = c(outcomes1 , hat$minimum)

hx=N^(-1/5)
ht=N^(-1/5)
hat=optimize(MHD,lower = 0,upper = 10)
outcomes2 = c(outcomes2 , hat$minimum)


fit = coxph(Surv(alltime,del)~Z)
outcomes3 = c(outcomes3,fit$coefficients[1])
} 
mean(outcomes1)- beta0   # cons_h0 vs est_h0  bandwith = n^(-1/6)
sum((outcomes1-mean(outcomes1))^2)/length(outcomes1)

mean(outcomes2)- beta0   # cons_h0 vs est_h0  bandwith = n^(-1/5)
sum((outcomes2-mean(outcomes2))^2)/length(outcomes2)

mean(outcomes3)- beta0   #PMLE
sum((outcomes3-mean(outcomes3))^2)/length(outcomes3)
proc.time()-ptm