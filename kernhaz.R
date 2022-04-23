rm(list=ls())
install.packages("kernhaz")
library(kernhaz)

beta0 = 2
N=300   
h=(3*N)^(-1/3)
b=1
library(survival)
known_h0 = function(t)0.1+0.2*t
plot(xx,known_h0(xx))

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
  
  rhat <-function(y){
    s=0
    for (i in 1:length(alltime)){
      s= s+wstar0[i]*kb(y-alltime[i])
    }
    s/sum(wstar0)
  }
  
  capK <- function(t){
    integrate(k,-Inf,t)$value
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
  
  
  rhat <-function(y){
    s=0
    for (i in 1:length(alltime)){
      s= s+wstar0[i]*kb(y-alltime[i])
    }
    s/sum(wstar0)
  }
  capK <- function(t){
    integrate(k,-Inf,t)$value
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



BLestimator <- function (Z,delta,betaZ) 
{ UZ<-unique(Z)      #Unique observed times;
N<-length(UZ)
UZ.order<-order(UZ)
UZ<-UZ[UZ.order]            #sort data;
KM <- rep(0,N)
BrelowH <- rep(0,N)
EBZ=exp(betaZ)

Y<-rep(0,N)
D<-rep(0,N)
h<-rep(0,N)
SEBZ<-rep(0,N)
D[1]<-sum(Z[delta==1]==UZ[1])
Y[1]<-sum(Z >= UZ[1])
SEBZ[1] <-sum((Z >= UZ[1])*EBZ)


BrelowH[1] <-D[1]/SEBZ[1]
h[1]=BrelowH[1]
for (i in 2: N){
  D[i]<-sum(Z[delta==1]==UZ[i])
  Y[i]<-sum(Z >= UZ[i])
  SEBZ[i] <-sum((Z >= UZ[i])*EBZ)
  BrelowH[i] <-BrelowH[i-1]+D[i]/SEBZ[i]
  h[i]<-BrelowH[i]-BrelowH[i-1]
}
data.frame(UZ,D,Y,BrelowH)
}



MHD2 <-function(b0){
  
  b0 = beta0
  zeroset = rep(0,length(alltime))
  for (i in 1:length(alltime)){
    zeroset[i]=condition_hazard(alltime,Z[i],Z)(alltime[i])/exp(b0*Z[i])
  }
  
  cons_h0 = function(t)  kernelreg(alltime,zeroset,t)
  est_h0 = condition_0hazard(alltime,Z)
  plot(xx,cons_h0(xx))
  
  distance = 0
  for (i in 1:length(alltime)){
    distance = distance + (cons_h0(alltime[i])^0.5-est_h0(alltime[i])^0.5)^2
    
  }
  distance
}



  newdata0 = inversion_gen(N/3,0)  # 0-8
  newdata1 = inversion_gen(N/3,1)  # 0-6
  newdata2 = inversion_gen(N/3,2)  # 0-4
  alltime = c(newdata0,newdata1,newdata2)
  
  Z= c(rep(0,N/3),rep(1,N/3),rep(2,N/3))
  
  
  Z = runif(N,0.5,5)
  alltime= c()
  for (i in 1:N){
    temp = inversion_gen(1,Z[i])
    alltime = c(alltime,temp)
  }
  alldata = cbind(alltime,Z)
  alldata = alldata[order(alldata[,1]),]
  alltime = alldata[,1]
  Z = alldata[,2]
  del = rep(1,N)

  selfcon_h0 = condition_0hazard(alltime,Z)
  xx= seq(0,max(alltime),0.1)
  yy= rep(0,length(xx))
  for (i in 1:length(xx)){
    yy[i]=selfcon_h0(xx[i])
  }
  plot(xx,yy)
  
  r = khazardcond(alltime,del,Z,h=c(N^(-1/4),1),t=alltime,x=Z)
  plot(r$time.points,r$hazard)
  condhazard = r$hazard
  plot(r$time.points,  condhazard[1,])
  rcon_h0 = function(t)kernelreg(r$time.points,r$hazard,t)
  plot(xx,rcon_h0(xx))  
  khazardcond(alltime,del,Z,h=c(5,1),t=alltime[i],x=Z[i])$hazard
  condition_hazard(alltime,Z[i],Z)(alltime[i])
  i=2
  