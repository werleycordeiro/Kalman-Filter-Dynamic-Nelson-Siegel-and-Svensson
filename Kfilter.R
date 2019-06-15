# Filter
Kfilter<-function(logLik,N,T,Y,Z,a.t,P.t,H,a.tt,P.tt,v2,v1,phi,mu,Q,prev,M,Yf,lik){
for (t in 1:T) 
{
  v <- (as.numeric(Y[t, ])) - Z %*% a.t[t, ] # prediciton error vector
  F <- Z %*% P.t[t, ,] %*% t(Z) + H # prediciton error variance matrix
  if(det(F)<=1e-30 || is.na(det(F)) || is.nan(det(F)) || is.infinite(det(F))){
    logLik<- -1000000000000000; break
  }else{
    F.inv  <- solve(F)
    # Log-likelihood
    logLik <- logLik - 0.5 * (log(det(F)) + t(v) %*% F.inv %*% v) # constructed via the prediction error decomposition
    }
    # Updating the state vector and its variance matrix
    a.tt[t, ]   <- a.t[t, ] +  P.t[t, , ] %*% t(Z) %*% F.inv %*% v
    P.tt[t, , ] <- P.t[t, , ] - P.t[t, , ] %*% t(Z) %*% F.inv %*% Z %*% P.t[t, , ]
    v1[t, ]	<- Z %*% a.tt[t, ]
    v2[t, ] <- (as.numeric(Y[t, ])) - Z %*% a.tt[t, ] # Filtered errors
    # Predicting the state vector and its variance matrix
    a.t[t + 1, ]  <- phi %*% a.tt[t, ] + (diag(N) - phi) %*% mu  
    P.t[t + 1, ,] <- phi %*% P.tt[t, ,] %*% t(phi) + Q
  }
  # Forecast
  if(prev) 
  {
    if(t > (T - 1))
      for(m in 1:M)
      {
        Yf[t + m,]<- Z %*% a.t[t + m, ] 
        a.tt[t + m, ]<- a.t[t + m, ]
        P.tt[t + m, , ]<- P.t[t + m, , ]
        a.t[t + m + 1, ]  <- phi %*% a.tt[t + m, ] + (diag(N) - phi) %*% mu  
        P.t[t + m + 1, ,] <- phi %*% P.tt[t + m, ,] %*% t(phi) + Q
      }
  }
  if(lik)
  {
    as.numeric(-logLik)
  }else{
    return(list(a.tt=a.tt,a.t=a.t,P.tt=P.tt,P.t=P.t,v2=v2,v1=v1,Yf=Yf)) #****
  }
}	
