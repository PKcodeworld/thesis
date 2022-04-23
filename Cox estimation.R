rm(list=ls())

par(mfrow=c(1,1))
beta0 = 2
N=600
h=(3*N)^(-1/6)
b=(3*N)^(-1/6)
library(survival)

known_h0 = function(t)0.1+0.2*t

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


#############################
MHD1 <-function(b0){

  zeroset = rep(0,length(alltime))
  for (i in 1:length(alltime)){
    zeroset[i]=condition_hazard(alltime,Z[i],Z)(alltime[i])/exp(b0*Z[i])
  }

  cons_h0 = function(t)  kernelreg(alltime,zeroset,t)
 
  #cons_h0 = function(t)  kernelreg(ddsub[,1],ddsub[,2],t)
 # xx= seq(0,max(alltime),0.1)
#  plot(xx,cons_h0(xx))
 # abline(lr$coefficients[1],lr$coefficients[2],col="red")
 # abline(0.1,0.2,col="black")

 # plot(alltime,zeroset)

 # abline(lr$coefficients[1],lr$coefficients[2])
 # abline(0.1,0.2)
  #ddd=data.frame(alltime,zeroset)
 # ddsub = subset(ddd,ddd$alltime<=4)
# library(ggplot2)
#  lr = lm(zeroset~alltime)
 # plot(alltime,zeroset)
 # po = ggplot(ddd,aes(x=alltime,y=zeroset))+geom_point(shape = 20)
 # po+geom_smooth(method=lm)
 # summary(l)
  
 # betaZ = beta0*Z
 # del = rep(1,length(alltime))
 # l = BLestimator(alltime,del,betaZ)
 # curv= function(t) kernelreg(l[,1],l[,4],t)
#  dx = 0.0001
#  haz = function(t){(curv(t+dx)-curv(t))/dx}
#  est_h0 = function(t) kernelreg(alltime,haz(alltime),t)
#  plot(xx,est_h0(xx))
#  abline(0.1,0.2)
 # plot(xx,curv(xx))
  
#  condi_h0 = condition_hazard(alltime,0,Z)
  #yy= rep(0,length(xx))
  #for (i in 1:length(xx)){
  #  yy[i]=condi_h0(xx[i])
  #}
  #plot(xx,yy)
  #abline(0.1,0.2)
  est_h0= function (t)0.1+0.2*t  
  distance = 0
  for (i in 1:length(alltime)){
    distance = distance + (cons_h0(alltime[i])^0.5-est_h0(alltime[i])^0.5)^2
    
  }
  distance
}


bb=seq(0,5,0.1)
ddis=rep(0,length(bb))
for (i in 1:length(bb)){
  ddis[i] = MHD1(bb[i])
}

plot(bb,ddis)


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

  newdata0 = inversion_gen(N/3,3)  # 0-8
  newdata1 = inversion_gen(N/3,1)  # 0-6
  newdata2 = inversion_gen(N/3,2)  # 0-4
  alltime = c(newdata0,newdata1,newdata2)
  
  Z= c(rep(3,N/3),rep(1,N/3),rep(2,N/3))
  
  
  Z = runif(N,0.5,3)
  alltime= c()
  for (i in 1:N){
    temp = inversion_gen(1,Z[i])
    alltime = c(alltime,temp)
  }
  alldata = cbind(alltime,Z)
  alldata = alldata[order(alldata[,1]),]
  alltime = alldata[,1]
  Z = alldata[,2]

  ptm<-proc.time()
  hat1 = optimize(MHD1,lower = 0,upper = 10)
proc.time()-  ptm
hat1$minimum

library(muhaz)
del = rep(1,N)    #censoring indicator


fit0<-muhaz(newdata0,del,n.est.grid = N)
fit1<-muhaz(newdata1,del,n.est.grid = N)
fit2<-muhaz(newdata2,del,n.est.grid = N)
haz0fun <-function(t) {kernelreg(fit0$est.grid,fit0$haz.est,t)}
haz1fun <-function(t) {kernelreg(fit1$est.grid,fit1$haz.est,t)}
haz2fun <-function(t) {kernelreg(fit2$est.grid,fit2$haz.est,t)}
alltime1 = c(fit0$est.grid,fit1$est.grid,fit2$est.grid)

all_0hazard1=c(fit0$haz.est/exp(0*beta0),fit1$haz.est/exp(1*beta0),fit2$haz.est/exp(2*beta0))
plot(alltime1,all_0hazard1)
bline  <-function(t) {kernelreg(alltime1,all_0hazard1,t)}

plot(xx,bline(xx))
plot(fit2$est.grid,fit2$haz.est)

