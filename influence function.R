############# initial setting ########
#h0(t)=0.1+0.2t
#Z = 0,1,2
rm(list=ls())
beta0 = 2


#____________generate data_____________________
N=50   # 50 - 100 - 300

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
    numer = numer+y[i]*K((t-x[i])/h)
    deno= deno+K((t-x[i])/h)
  }
  return (numer/deno)
  
}
#########################


Hellinger_distance <-function(beta){
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

influenceH=c()
influenceP= c()
for(i in 1:100){
newdata0 = inversion_gen(N,0)  # 0-8
newdata1 = inversion_gen(N,1)  # 0-6
newdata2 = inversion_gen(N,2)  # 0-4
newdata2 = sort(newdata2)

# setting outliners
del = rep(1,N)    #censoring indicator
fit0<-muhaz(newdata0,del,n.est.grid = N)
fit1<-muhaz(newdata1,del,n.est.grid = N)
fit2<-muhaz(newdata2,del,n.est.grid = N)

hat = optimize(Hellinger_distance,lower = 0,upper = 5)
H_origin_beta = hat$minimum

alltime = c(newdata0,newdata1,newdata2)
alldel = rep(1,(3*N))
g = c(rep(0,N),rep(1,N),rep(2,N))
coxfit = coxph(Surv(alltime,alldel)~g)
P_origin_beta=coxfit$coefficients[1]


zz=seq(0,20,0.1)
influH = c()
influP = c()
for (i in 1:length(zz)){
  newdata2[1]=rnorm(1,zz[i],0.1)
  
  fit2<-muhaz(newdata2,del,n.est.grid = N)
  hat2=optimize(Hellinger_distance,lower = 0,upper = 5)
  H_contaminated_beta=hat2$minimum
  
  alltime = c(newdata0,newdata1,newdata2)
  alldel = rep(1,3*N)
  g = c(rep(0,N),rep(1,N),rep(2,N))
  coxfit = coxph(Surv(alltime,alldel)~g)
  P_contaminated_beta=coxfit$coefficients[1]

  influH=c(influH,N*(H_contaminated_beta-H_origin_beta))
  influP=c(influP,N*(P_contaminated_beta-P_origin_beta))
}
influenceH=c(influenceH,influH)
influenceP=c(influenceP,influP)

}
H=matrix(influenceH,ncol=length(zz),byrow=TRUE)
P=matrix(influenceP,ncol=length(zz),byrow=TRUE)


par(mfrow=c(1,2))
Hsum=colSums(H)
Hmean=Hsum/100
plot(zz,Hmean)
Psum=colSums(P)
Pmean=Psum/100
plot(zz,Pmean)

