rm(list=ls())

beta0 = 0.5
N=50   
h=(3*N)^(-1/6)
b=(3*N)^(-1/6)

inversion_gen<-function(n,Z){
  temp = runif(n)
  x = (-1+sqrt(1-40*log(1-temp)/exp(beta0*Z)))/2
}

newdata0 = inversion_gen(N,0)  # 0-8
newdata1 = inversion_gen(N,1)  # 0-6
newdata2 = inversion_gen(N,2)  # 0-4
alltime = c(newdata0,newdata1,newdata2)

Z= c(rep(0,N),rep(1,N),rep(2,N))

condition_hazard <- function(alltime,Z){

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

hazestimate = condition_hazard(alltime,Z)


xx= seq(0,4,0.01)
yy= rep(0,length(xx))
for (i in 1:length(xx)){
  yy[i]=hazestimate(xx[i])  
}


plot(xx,yy)

