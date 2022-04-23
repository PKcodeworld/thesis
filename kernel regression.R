############# initial setting ########
beta0 = 0.5
#h0(t)=0.1+0.2t
#Z = 0,1,2
rm(list=ls())
#############  test area #############
p <- function(x){
  (0.1+0.2*x)*exp(-(0.1*x+0.1*x^2))
}
p2<-function(x){
  exp(0.5)*(0.1+0.2*x)*exp(-exp(0.5)*(0.1*x+0.1*x^2))
}

hist(newdata1)

x = seq(0,8,0.01)

par(mfrow=c(1,1))
plot(x,p2(x))
p2(0)

xx = seq(0,8,0.01)
plot(xx,p(xx))

plot(fit$time,fit$hazard)
max(fit$time)

rm(list=ls())
#____________generate data_____________________
N=500   # 50 - 100 - 300
# censoring rate 10%-30%
# repeat 500 times calculate MSE/Bias
inversion_gen<-function(n,Z){
  temp = runif(n)
  x = (-1+sqrt(1-40*log(1-temp)/exp(beta0*Z)))/2
}

newdata0 = inversion_gen(N,0)  # 0-8
newdata1 = inversion_gen(N,1)  # 0-6
newdata2 = inversion_gen(N,2)  # 0-4

#______hazard estimation based on interval-Poisson assumption_______

install.packages("bshazard")
library(bshazard)

del = rep(1,N)    #censoring indicator
fit0<-bshazard(Surv(newdata0,del)~1)
fit1<-bshazard(Surv(newdata1,del)~1)
fit2<-bshazard(Surv(newdata2,del)~1)

plot(fit2$time,fit2$hazard)
xx=seq(0,10,0.01)
haz = rep(0,length(xx))
for (i in 1:length(xx)){
  haz[i] = haz2fun(xx[i])
}
plot(xx,haz)

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

Hellinger_distance <-function(beta){
  b = beta
  alltime = c(fit0$time,fit1$time,fit2$time)
  all_0hazard=c(fit0$hazard/exp(0*b),fit1$hazard/exp(1*b),fit2$hazard/exp(2*b))
  
  hazfun <- function(t){
    kernelreg(alltime,all_0hazard,t)
  }
  
  distance0 <- function(x){
    (sqrt(haz0fun(x))-sqrt(hazfun(x)))^2
  }

  distance1 <- function(x){
    (sqrt(haz1fun(x))-sqrt(hazfun(x)*exp(b*1)))^2
  }

  
  distance2 <- function(x){
    (sqrt(haz2fun(x))-sqrt(hazfun(x)*exp(b*2)))^2
  }
  
 # x=seq(0,10,0.1)
 # y=rep(0,length(x))
 # for (i in 1:length(x)){
  #  y[i]=distance0(x[i])
 # }
 # plot(x,y)
  low0=min(newdata0)
  low1=min(newdata1)
  low2=min(newdata2)
  
  up0 = max(newdata0)
  up1 = max(newdata1)
  up2 = max(newdata2)
  integrate(distance0,low0,up0)$value+integrate(distance1,low1,up1)$value+integrate(distance2,low2,up2)$value
}


b = seq(0,0.7,0.01)
dis = rep(0,length(b))
for (i in 1:length(b)){
  dis[i] = Hellinger_distance(b[i])
}
plot(b,dis)

hat = optimize(Hellinger_distance,lower = 0,upper = 0.7)
hat
## ____ plug in the bhat to estimate lambda0 (baseline hazard)____
alltime = c(fit0$time,fit1$time,fit2$time)
alldel = rep(1,3*N)
bhat = hat$minimum
all_0hazard=c(fit0$hazard/exp(0*bhat),fit1$hazard/exp(1*bhat),fit2$hazard/exp(2*bhat))
lambda0 = lm(all_0hazard~alltime)
summary(lambda0)


##spline regression__________________

library(splines)
knots = quantile(fit0$time,c(0.25,0.5,0.75))
model0_gam = lm(fit0$hazard~bs(fit0$time,knots=knots))
summary(model0_gam)
plot(model0_gam)

#________________guassian kernel/changeable_____________
NAestimator <- function (Z,delta) 
{ 
 
  UZ<-unique(Z)      #Unique observed times;
N<-length(UZ)
UZ.order<-order(UZ)
UZ<-UZ[UZ.order]            #sort data;
N.A <- rep(0,N)             #cannot use NA (reserved by R) for N.A;
Y<-rep(0,N)
D<-rep(0,N)
D[1]<-sum(Z[delta==1]==UZ[1])
Y[1]<-sum(Z >= UZ[1])
N.A[1] <- D[1]/Y[1]            #this is for right continuous value
for (i in 2: N){
  D[i]<-sum(Z[delta==1]==UZ[i])
  Y[i]<-sum(Z >= UZ[i])
  N.A[i] <- N.A[i-1]+D[i]/Y[i]
}

# Calculate variance and sdandard error of N-A estimator;
h = rep(0,N)
h[1]=N.A[1]
for (i in 2: N){
  h[i]=N.A[i]-N.A[i-1]
}


Full<-data.frame(UZ,D,N.A,h)
Reduced<-subset(Full, (Full$D>0))
Reduced
}
fit0= NAestimator(newdata0,del)
fit1= NAestimator(newdata1,del)
fit2= NAestimator(newdata2,del)
colnames(fit0)<-c("time","death","hazard","cumulH")
colnames(fit1)<-c("time","death","hazard","cumulH")
colnames(fit2)<-c("time","death","hazard","cumulH")
fit0
#########################################

alltime = c(fit0$time,fit1$time,fit2$time)
alldel = rep(1,3*N)
g = c(rep(0,N),rep(1,N),rep(2,N))
library(survival)
coxfit = coxph(Surv(alltime,alldel)~g)
coxfit$coefficients[1]
