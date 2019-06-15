# Filter
Kfilter<-function(N,T,Y,Z,a.t,P.t,H,a.tt,P.tt,v2,phi,mu,Q,prev,M,Yf){
for (t in 1:T) 
{
  v <- (as.numeric(Y[t, ])) - pars$Z %*% a.t[t, ] # prediciton error vector
  F <- pars$Z %*% P.t[t, ,] %*% t(pars$Z) + H # prediciton error variance matrix
  if(det(F)<=1e-30 || is.na(det(F)) || is.nan(det(F)) || is.infinite(det(F))){
    logLik<- -1000000000000000; break
  }else{
    F.inv  <- solve(F)
    # Log-likelihood
    logLik <- logLik - 0.5 * (log(det(F)) + t(v) %*% F.inv %*% v) # constructed via the prediction error decomposition
    # Updating the state vector and its variance matrix
    a.tt[t, ]   <- a.t[t, ] +  P.t[t, , ] %*% t(pars$Z) %*% F.inv %*% v
    P.tt[t, , ] <- P.t[t, , ] - P.t[t, , ] %*% t(pars$Z) %*% F.inv %*% pars$Z %*% P.t[t, , ]
    v2[t, ] <- (as.numeric(Y[t, ])) - pars$Z %*% a.tt[t, ] # Filtered errors
    # Predicting the state vector and its variance matrix
    a.t[t + 1, ]  <- pars$phi %*% a.tt[t, ] + (diag(N) - pars$phi) %*% pars$mu  
    P.t[t + 1, ,] <- pars$phi %*% P.tt[t, ,] %*% t(pars$phi) + Q
  }
  # Forecast
  if(prev) 
  {
    if(t > (T - 1))
      for(m in 1:M)
      {
        Yf[t + m,]<- pars$Z %*% a.t[t + m, ] 
        a.tt[t + m, ]<- a.t[t + m, ]
        P.tt[t + m, , ]<- P.t[t + m, , ]
        a.t[t + m + 1, ]  <- pars$phi %*% a.tt[t + m, ] + (diag(N) - pars$phi) %*% pars$mu  
        P.t[t + m + 1, ,] <- pars$phi %*% P.tt[t + m, ,] %*% t(pars$phi) + Q
      }
  }
}	
}