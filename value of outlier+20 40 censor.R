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

MHD40 = function(b){
  
  zeroset40= rep(0,length(alltime40))
  for (i in 1:length(alltime40)){
    zeroset40[i]=condional_haz(Z40,alltime40,del40,alltime40[i],Z40[i])/exp(b*Z40[i])
  }
  
  cons_h0 = function(t){kernelreg(alltime40,zeroset40,t)}
  
  #est_h0= function(t){condional_haz(Z,alltime,del,t,0)}
  est_h0 = function(t) 0.1+0.2*t
  distance= 0
  for (i in 1:length(alltime40)){
    distance = distance+(cons_h0(alltime40[i])^0.5-est_h0(alltime40[i])^0.5)^2
  }
  distance
}

H=c()
P=c()
H40=c()
P40=c()
numcol = 50
iterate_time = 10
ptm<-proc.time()
for (i in 1:iterate_time){
  Z = runif(N,0.2,2)
  alltime= c()
  for (i in 1:N){
    temp = inversion_gen(1,Z[i])
    alltime = c(alltime,temp)
  }
  del = rep(1,N)
  
  Z40=Z
  alltime40=alltime
  del40=del
  # 20% censor
  for (i in 1:N){
    temp = runif(1,0,8.4)
    if (temp<alltime[i]){
      alltime[i]=temp
      del[i]=0
    }
  }
  # 40% censor
  for (i in 1:N){
    temp = runif(1,0,4.2)
    if (temp<alltime40[i]){
      alltime40[i]=temp
      del40[i]=0
    }
  }
  
  hx=N^(-1/5)
  ht=N^(-1/5)

  influH = c()
  influH40 = c()
  influP=c()
  influP40=c()
  
  index20 = sample(which(del==1),1)
  index40 = sample(which(del40==1),1)
  
  xx=seq(0.1,20,4)
  for (x in xx){
    
    Z[index20]=x
    Z40[index40]=x
    
    hat=optimize(MHD,lower = 0,upper = 10)
    H_contaminated_beta=hat$minimum
    
    coxfit = coxph(Surv(alltime,del)~Z)
    P_contaminated_beta=coxfit$coefficients[1]
    
    influH=c(influH,H_contaminated_beta)
    influP=c(influP,P_contaminated_beta)
    
    hat40=optimize(MHD40,lower = 0,upper = 10)
    H_contaminated_beta40=hat40$minimum
    
    coxfit40 = coxph(Surv(alltime40,del40)~Z40)
    P_contaminated_beta40=coxfit40$coefficients[1]
    
    influH40=c(influH40,H_contaminated_beta40)
    influP40=c(influP40,P_contaminated_beta40)
  }
  H=c(H,influH)
  P=c(P,influP)
  H40=c(H40,influH40)
  P40=c(P40,influP40)
  
} 


proc.time()-ptm

H = matrix(H,ncol=length(xx),byrow=TRUE)
P = matrix(P,ncol=length(xx),byrow=TRUE)
H40 = matrix(H40,ncol=length(xx),byrow=TRUE)
P40 = matrix(P40,ncol=length(xx),byrow=TRUE)

Hsum = colSums(H)
Hmean = Hsum/iterate_time
Hbias = (Hmean-beta0)*N

Psum = colSums(P)
Pmean = Psum/iterate_time
Pbias = (Pmean - beta0)*N

H40sum = colSums(H40)
H40mean = H40sum/iterate_time
H40bias = (H40mean-beta0)*N

P40sum = colSums(P40)
P40mean = P40sum/iterate_time
P40bias = (P40mean - beta0)*N


plot(xx,Hbias,col=1,type="b",pch=21,ylim=c(-1*N,0.5*N),xlab = "value of outliers",ylab ="bias*N")

lines(xx,Pbias,col=2,type="b",pch=22)
lines(xx,H40bias,col=3,type="b",,pch=23)
lines(xx,P40bias,col=4,type="b",,pch=24)
legend("topright",x.intersp=1,y.intersp=0.2,pch=21:24,text.width=5,legend=c("Hbias(20%censor)","Pbias(20%censor)","Hbias(40%censor)","Pbias(40%censor)"),col=c(1,2,3,4),bty="n")

