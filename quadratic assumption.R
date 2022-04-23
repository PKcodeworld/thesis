############# initial setting ########
beta0 = 0.5
#h0(t)=0.1+0.2t+0.3t^2
#Z = 0,1,2

#############  test area #############
f0 <- function(x){
  (0.1+0.2*x+0.3*x^2)*exp(-(0.1*x+0.1*x^2+0.1*x^3))
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
par(mfrow=c(1,2))

xx = seq(0,8,0.01)
plot(xx,f0(xx))
plot(xx,f1(xx))
plot(xx,f2(xx))

plot(fit$time,fit$hazard)
max(fit$time)

#____________generate data_____________________
N=1000

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

install.packages("bshazard")
library(bshazard)

del = rep(1,N)    #censoring indicator
fit0<-bshazard(Surv(newdata0,del)~1)
fit1<-bshazard(Surv(newdata1,del)~1)
fit2<-bshazard(Surv(newdata2,del)~1)

plot(fit0$time,fit0$hazard)

model0 = lm(fit0$hazard~fit0$time+I((fit0$time)^2))
model1 = lm(fit1$hazard~fit1$time+I((fit1$time)^2))
model2 = lm(fit2$hazard~fit2$time+I((fit2$time)^2))

haz0fun <- function(x){
  model0$coefficients[1]+model0$coefficients[2]*x+model0$coefficients[3]*x^2
}

haz1fun <- function(x){
  model1$coefficients[1]+model1$coefficients[2]*x+model1$coefficients[3]*x^2
}
haz2fun <- function(x){
  model2$coefficients[1]+model2$coefficients[2]*x+model2$coefficients[3]*x^2
}



Hellinger_distance <-function(beta){
  b = beta
  
  alltime = c(fit0$time,fit1$time,fit2$time)
  all_0hazard=c(fit0$hazard/exp(0*b),fit1$hazard/exp(1*b),fit2$hazard/exp(2*b))
  
  
  allmodel = lm(all_0hazard~alltime+I(alltime^2))
  #kernel regression(non-parametric,cr rao)
  
  hazfun <- function(x){
    allmodel$coefficients[1]+allmodel$coefficients[2]*x+allmodel$coefficients[3]*x^2
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

  
  integrate(distance0,0,100)$value+integrate(distance1,0,100)$value+integrate(distance2,0,100)$value
}


t = seq(0,0.7,0.001)
dis = rep(0,length(t))
for (i in 1:length(t)){
  dis[i] = Hellinger_distance(t[i])
}
plot(t,dis)

hat = optimize(Hellinger_distance,lower = 0,upper = 0.7)
hat
## ____ plug in the bhat to estimate lambda0 (baseline hazard)____
alltime = c(fit0$time,fit1$time,fit2$time)
bhat = hat$minimum
all_0hazard=c(fit0$hazard/exp(0*bhat),fit1$hazard/exp(1*bhat),fit2$hazard/exp(2*bhat))
lambda0 = lm(all_0hazard~alltime+I(alltime^2))
summary(lambda0)


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
