rm(list=ls())
library(KMsurv)
library(survival)
data(hodg)
burn = hodg

#gtype: Graft type (1=allogenic, 2=autologous) 
#dtype: Disease type (1=Non Hodgkin lymphoma, 2=Hodgkins disease) 
#time: Time to death or relapse, days 
#delta: Death/relapse indicator (0=alive, 1=dead) 
#rescale the covariates
#burn$time = (burn$time-min(burn$time))/(max(burn$time)-min(burn$time))
burn$time = burn$time/30
burn$gtype = burn$gtype-1  #gtype: Graft type (0=allogenic, 1=autologous) 
burn$dtype = burn$dtype-1 #dtype: Disease type (0=Non Hodgkin lymphoma, 1=Hodgkins disease) 
data0 = subset(burn,gtype ==0&dtype==0)
data1 = subset(burn,gtype ==0&dtype==1)
data2 = subset(burn,gtype ==1&dtype==0)
data3 = subset(burn,gtype ==1&dtype==1)


#kernel regression  ##############################################
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

######################################################

coxfit = coxph(Surv(time,delta)~gtype+dtype,data=burn)
summary(coxfit)
coxfit$coefficients
#  gtype   dtype 
# -0.2599  0.2742 



N = nrow(burn)
hazard = rep(0,N)
zeroset = rep(0,N)
for (i in 1:N){
  hazard[i] = condional_haz(burn$gtype,burn$dtype,burn$time,burn$delta,burn$time[i],burn$gtype[i],burn$dtype[i])
}

baseline = function(t){condional_haz(burn$gtype,burn$dtype,burn$time,burn$delta,t,0,0)}

Hellinger_distance <-function(beta){
  b0 = beta[1]
  b1 = beta[2]
  
  for (i in 1:N){
    zeroset[i]=hazard[i]/exp(b0*burn$gtype[i]+b1*burn$dtype[i])
  }
  
  sum = 0
  for (i in 1:N){
    sum = sum + (baseline(burn$time[i])^(1/2)-zeroset[i]^(1/2))^2
  }
  sum
}
es=optim(c(0,0),Hellinger_distance)
es$par    
