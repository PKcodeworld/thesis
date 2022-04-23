rm(list=ls())

beta0 = 2
N=900   
h=(3*N)^(-1/6)
b=(3*N)^(-1/6)
library(survival)

inversion_gen<-function(n,Z){
  temp = runif(n)
  x = (-1+sqrt(1-40*log(1-temp)/exp(beta0*Z)))/2
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
condition_0hazard <- function(alltime,Z){
    k <- function(x)  {1/sqrt(2*pi) * exp(-(x^2)/2)}
  
  w0 = rep(0,length(Z))
  for (i in 1:length(Z)){
    w0[i] = k(-Z[i]/h)/h
  }
  
  wstar0 = rep(0,length(Z))
  for (i in 1:length(Z)){
    wstar0[i] = w0[i]*(sum(w0*Z^2)-Z[i]*sum(w0*Z))
  }
  
  kb <-function(x){k(x/b)/b}
  
  rhat<-function(y){
    sum(wstar0*kb(y-alltime))/sum(wstar0)
  }
  
  capK <- function(y){
    integrate(k,-Inf,y)$value
  }
  
  Kb <-function(y){capK(y/b)}
  
  Hhat <-function(y){
    s= 0
    for (i in 1:length(wstar0)){
      s=s+wstar0[i]*Kb(alltime[i]-y)
    }
    s/sum(wstar0)
  }
  
  hazhat <- function(y){
    rhat(y)/Hhat(y) 
  }
  hazhat
}


condition_hazard <- function(alltime,z0,Z){
  k <- function(x)  {1/sqrt(2*pi) * exp(-(x^2)/2)}
  
  w0 = rep(0,length(Z))
  for (i in 1:length(Z)){
    w0[i] = k((z0-Z[i])/h)/h
  }
  
  wstar0 = rep(0,length(Z))
  for (i in 1:length(Z)){
    wstar0[i] = w0[i]*(sum(w0*(z0-Z)^2)-(z0-Z[i])*sum(w0*(z0-Z)))
  }
  
  kb <-function(x){k(x/b)/b}
  
  rhat<-function(y){
    sum(wstar0*kb(y-alltime))/sum(wstar0)
  }
  
  capK <- function(y){
    integrate(k,-Inf,y)$value
  }
  
  Kb <-function(y){capK(y/b)}
  
  Hhat <-function(y){
    s= 0
    for (i in 1:length(wstar0)){
      s=s+wstar0[i]*Kb(alltime[i]-y)
    }
    s/sum(wstar0)
  }
  
  hazhat <- function(y){
    rhat(y)/Hhat(y) 
  }
  hazhat
}

############################
#h2 = condition_hazard(alltime,2,Z)

#hazestimate = condition_0hazard(alltime,Z)
#xx= seq(0,2,0.01)
#yy= rep(0,length(xx))
#for (i in 1:length(xx)){
#  yy[i]=h2(xx[i])  
#}

#plot(xx,yy)
#############################
MHD <-function(b){
zeroset = rep(0,length(alltime))
for (i in 1:length(alltime)){
  zeroset[i]=condition_hazard(alltime,Z[i],Z)(alltime[i])/exp(b*Z[i])
}

cons_h0 = function(t)  kernelreg(alltime,zeroset,t)
#est_h0 = function(t)0.1+0.2*t
est_h0 = condition_0hazard(alltime,Z)
distance = 0
for (i in 1:length(alltime)){
distance = distance + (cons_h0(alltime[i])^0.5-est_h0(alltime[i])^0.5)^2

}
distance}

outcomes1 = c()
outcomes2 = c()
ptm<-proc.time()
for (i in 1:100)
{
#newdata0 = inversion_gen(N,0)  # 0-8
#newdata1 = inversion_gen(N,1)  # 0-6
#newdata2 = inversion_gen(N,2)  # 0-4
#alltime = c(newdata0,newdata1,newdata2)

#Z= c(rep(0,N),rep(1,N),rep(2,N))

  
  Z = runif(N,0,2)
alltime= c()
for (i in 1:N){
  temp = inversion_gen(1,Z[i])
  alltime = c(alltime,temp)
}

hat = optimize(MHD,lower = 0,upper = 10)


outcomes1 = c(outcomes1 , hat$minimum)


del = rep(1,N)
fit = coxph(Surv(alltime,del)~Z)
outcomes2 = c(outcomes2,fit$coefficients[1])
}

mean(outcomes1)- beta0   #0/1/2 cons_h0 vs est_h0
sum((outcomes1-mean(outcomes1))^2)/length(outcomes1)

mean(outcomes2)- beta0   #PMLE
sum((outcomes2-mean(outcomes2))^2)/length(outcomes2)
proc.time()-ptm