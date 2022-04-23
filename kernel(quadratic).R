############# initial setting ########
beta0 = 0.5
#h0(t)=0.1+0.2t+0.3t^2
#Z = 0,1,2
rm(list=ls())
#############  test area #############
f0 <- function(x){
  (0.1+0.2*x+0.3*x^2)*exp(-(0.1*x+0.1*x^2+0.1*x^3))
}
h0 <-function(t){
  0.1+0.2*t+0.3*t^2
}

f1<-function(x){
  (0.1+0.2*x+0.3*x^2)*exp(-exp(beta0)*(0.1*x+0.1*x^2+0.1*x^3))
}
f2<-function(x){
  (0.1+0.2*x+0.3*x^2)*exp(-exp(2*beta0)*(0.1*x+0.1*x^2+0.1*x^3))
}


hist(newdata0,freq = FALSE)
hist(newdata1,freq = FALSE)
hist(newdata2,freq = FALSE)
par(mfrow=c(1,1))

xx = seq(0,4,0.01)
plot(xx,f0(xx))
plot(xx,f1(xx))
plot(xx,f2(xx))

plot(fit$time,fit$hazard)
max(fit$time)

#____________generate data_____________________
N=100

inversion_gen<-function(n,Z){
  temp = runif(n)
  runif(1)
  a = -log(1-temp)/exp(beta0*Z)
  x = (3*sqrt(3)*sqrt(2700*a^2+140*a+3)+270*a+7)^(1/3)/(3*2^(1/3))-2*2^(1/3)/(3*(3*sqrt(3)*sqrt(2700*a^2+140*a+3)+270*a+7)^(1/3))-1/3
}


newdata0 = inversion_gen(N,0)
newdata1 = inversion_gen(N,1)
newdata2 = inversion_gen(N,2)

#______hazard estimation based on interval-Poisson assumption_______

library(bshazard)

del = rep(1,N)    #censoring indicator
fit0<-bshazard(Surv(newdata0,del)~1)
fit1<-bshazard(Surv(newdata1,del)~1)
fit2<-bshazard(Surv(newdata2,del)~1)

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
    (sqrt(haz1fun(x)/exp(b*1))-sqrt(hazfun(x)))^2
  }
  
  
  distance2 <- function(x){
    (sqrt(haz2fun(x)/exp(b*2))-sqrt(hazfun(x)))^2
  }
  
  # x=seq(0,10,0.1)
  # y=rep(0,length(x))
  # for (i in 1:length(x)){
  #  y[i]=distance0(x[i])
  # }
  # plot(x,y)
  
  up0 = max(fit0$time)
  low0 = min(fit0$time)
  up1 = max(fit1$time)
  low1 = min(fit1$time)
  up2 = max(fit2$time)
  low2 = min(fit2$time)
  integrate(distance0,low0,up0)$value+integrate(distance1,low1,up1)$value+integrate(distance2,low2,up2)$value
}

b = seq(0,0.7,0.01)
dis = rep(0,length(b))
for (i in 1:length(b)){
  dis[i] = Hellinger_distance(b[i])
}
plot(b,dis)

hat = optimize(Hellinger_distance,lower = 0,upper = 2)
hat$minimum-beta0
## ____ plug in the bhat to estimate lambda0 (baseline hazard)____
alltime = c(newdata0,newdata1,newdata2)
alldel = rep(1,3*N)
g = c(rep(0,N),rep(1,N),rep(2,N))
library(survival)
coxfit = coxph(Surv(alltime,alldel)~g)

hat$minimum-beta0
coxfit$coefficients[1]-beta0

##spline regression__________________

library(splines)
knots = quantile(fit0$time,c(0.25,0.5,0.75))
model0_gam = lm(fit0$hazard~bs(fit0$time,knots=knots))
summary(model0_gam)
plot(model0_gam)

#________________guassian kernel/changeable_____________

kernelplot<- function(x){
  
  n <- length(x)
  gauss <- function(x) 1/sqrt(2*pi) * exp(-(x^2)/2)
  xgrid <- seq(from = min(x) - 1, to = max(x) + 1, by = 0.01)
  h= 4*sd(x)/3*n^(-0.2)  #Silverman's Rule of Thumb
  
  bumps <- sapply(x, function(a) gauss((xgrid - a)/h)/(n *h))
  #(grid*sample size)
  
  plot(xgrid, rowSums(bumps), ylab = expression(hat(f)(x)),type = "l", xlab = "x", lwd = 2)
  out = data.frame(xgrid,rowSums(bumps))
  out
}


kernelfun<-function(x,t){
  n <- length(x)
  gauss <- function(x) 1/sqrt(2*pi) * exp(-(x^2)/2)
  h= 4*sd(x)/3*n^(-0.2)  #Silverman's Rule of Thumb
  return (sum(gauss((t-x)/h)/(n *h)))
  
}


#____________construct v(b)_______________

v<-function(x,b,Z){
  x/exp(b*Z)
}

zero1 = v(fit1$hazard,0.5,1)
kernelplot(zero1)
kernelfun(zero1,5)
plot(fit1$time,v(fit1$hazard,0.5,1))
