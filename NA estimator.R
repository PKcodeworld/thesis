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


