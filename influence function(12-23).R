rm(list=ls())

beta0 = 1
N=100  
library(survival)


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


influenceH=c()
influenceP= c()
ptm<-proc.time()
for (i in 1:100){
  Z = runif(N,0.2,1)
  alltime= c()
  for (i in 1:N){
    temp = inversion_gen(1,Z[i])
    alltime = c(alltime,temp)
  }
  del = rep(1,N)
  
  hx=N^(-1/5)
  ht=N^(-1/5)
  hat=optimize(MHD,lower = 0,upper = 10)
  H_origin_beta = hat$minimum
  
  
  coxfit = coxph(Surv(alltime,del)~Z)
  P_origin_beta=coxfit$coefficients[1]
  
  index = sample(1:N,1)

  
  zz=seq(0,100,1)
  influH = c()
  influP = c()
  
  for (i in 1:length(zz)){
    Z[index]=zz[i]
    
    hat2=optimize(MHD,lower = 0,upper = 10)
    H_contaminated_beta=hat2$minimum
    
    coxfit = coxph(Surv(alltime,del)~Z)
    P_contaminated_beta=coxfit$coefficients[1]
    
    influH=c(influH,N*(H_contaminated_beta-H_origin_beta))
    influP=c(influP,N*(P_contaminated_beta-P_origin_beta))
  }
  influenceH=c(influenceH,influH)
  influenceP=c(influenceP,influP)
 
} 

proc.time()-ptm

H=matrix(influenceH,ncol=length(zz),byrow=TRUE)
P=matrix(influenceP,ncol=length(zz),byrow=TRUE)


par(mfrow=c(1,2))
Hsum=colSums(H)
Hmean=Hsum/N
plot(zz,Hmean)
Psum=colSums(P)
Pmean=Psum/N
plot(zz,Pmean)
