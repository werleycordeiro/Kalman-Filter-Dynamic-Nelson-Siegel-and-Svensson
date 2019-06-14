# author: Werley Cordeiro
# werleycordeiro@gmail.com

# Data
#setwd("C:\\...")
#data <- read.csv("bonds.csv",header = TRUE, sep = ";")
data<-read.csv("https://www.dropbox.com/s/inpnlugzkddp42q/bonds.csv?dl=1",header = TRUE, sep = ";")#***
require(xts)
data<- as.matrix(data[,which(names(data)=="M3"):which(names(data)=="M120")])
datas<- seq(as.Date("1972/1/1"), by = "month", length.out = 348)
data<- xts(data, order.by = datas)
data1<- xts(data, order.by = datas) #***


# Initial Parameters (See Dynamic-Nelson-Siegel to initial Parameters)  
# Lambda(1x1), H diag matrix (17x17), phi matrix (3x3), mu vector (3x1), Q matrix (3x3). Total: 36 parameters

para<-c(0.0609,

0.14170940,0.07289485,0.11492339,0.11120008,0.09055795,0.07672075,0.07222108,0.07076431,0.07012891,0.07267366,0.10624206,0.09029621,0.10374527,0.09801215,0.09122014,0.11794190,0.13354418,

 0.99010443,0.02496842,-0.002294319,
-0.02812401,0.94256154, 0.028699387,
 0.05178493,0.01247332, 0.788078795, 

8.345444,-1.572442,0.2029919,  

 0.3408764,
-0.07882772,0.62661018,
-0.21351036,-0.00425989,1.08802059)

# Kalman Filter

prev<- FALSE # TRUE to Forecast.
ahead<-12 # X-step ahead forecast
lik <- TRUE # TRUE to return the value of the loglikelihood function. FALSE to return parameters.

kalman <- function(para,Y,lik,prev,ahead) {#***
  l   <- para[1]
  Nelson.Siegel.factor.loadings <- function(l)
  {
    m <- c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
    column1 <- rep.int(1,length(m))
    column2 <- (1 - exp(-l * m))/(l * m)
    column3 <- column2 - exp(-l * m)
    
    lambmat <- cbind(column1,column2,column3)
    
    lambmat
  }  
  
  M<-ahead #***
  
	  if(prev){#*** Forecast
		  T <- nrow(Y)
		  Yf<-Y
		  Yf[(T-M+1):T,]<-NA
		  Y<-Y[1:(T-M),]
		  T <- nrow(Y)
		  }else{
		  T <- nrow(Y)}#***
  pars<-list()
  W <- ncol(Y)
  N <- 3
  
  pars$mu	<- matrix(NA,N,1) # Mean vector
  pars$phi	<- diag(N) # Vector Autoregressive coeffient matrix VAR(1)	
  pars$H	<- diag(ncol(Y)) # Variance matrix of residuals
  pars$Q	<- diag(N) # Transition covariance matrix of residuals
  pars$Z	<- Nelson.Siegel.factor.loadings(l) # Loading matrix
  
  # Variance matrix of residuals
  for(i in 1:17){
  pars$H[i,i]<-para[1 + i]
  }

  H <- pars$H^2
  
  # Vector autoregressive coeffient matrix: VAR(1)
  pars$phi[1,1] <- para[19]
  pars$phi[1,2] <- para[20]
  pars$phi[1,3] <- para[21]
  pars$phi[2,1] <- para[22]
  pars$phi[2,2] <- para[23]
  pars$phi[2,3] <- para[24]
  pars$phi[3,1] <- para[25]
  pars$phi[3,2] <- para[26]
  pars$phi[3,3] <- para[27]
  
  # Mean vector
  pars$mu[1]<-para[28]
  pars$mu[2]<-para[29]
  pars$mu[3]<-para[30]
  
  # Transition covariance matrix of residuals
  pars$Q[1,1] <- para[31]
  pars$Q[2,1] <- para[32]
  pars$Q[2,2] <- para[33]
  pars$Q[3,1] <- para[34]
  pars$Q[3,2] <- para[35]
  pars$Q[3,3] <- para[36]
  
  Q <- pars$Q %*% t(pars$Q) 
  
  v2   <- matrix(NA,T,W) # Filtered errors: are defined as the difference between the observed yield curve and its filtered estimate from KF
  
  # Forecast 	
  if(prev){#*** 
	  a.tt <- matrix(NA, (T+M), N)
	  a.t  <- matrix(NA, (T+M+1), N) # if prev=TRUE, always will be dim(a.t)[1]=348
	  P.tt <- array(NA, c((T+M), N, N))
	  P.t  <- array(NA, c((T+M+1), N, N))
  }else{
	  a.tt <- matrix(NA, T, N)
	  a.t  <- matrix(NA, (T+1), N)
	  P.tt <- array(NA, c(T, N, N))
	  P.t  <- array(NA, c((T+1), N, N))
  }#***
  
  # Start state vector and variance matrix
  a.t[1, ]  <- pars$mu # Start state vector: pars$at0
  P.t[1, ,] <- matrix(solve(diag(N^2) - kronecker(pars$phi,pars$phi)) %*% matrix(Q,(N^2),1),N,N)  # Start variance matrix: Lyapunov equation. pars$Pt0 # 
  
  # Initial loglikelihood	
  logLik <- - 0.5 * T * ncol(Y) * log(2 * pi)
  
  # Filter
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
  if(lik)
  {
    as.numeric(-logLik)
  }else
  {
    return(list(a.tt=a.tt,a.t=a.t,P.tt=P.tt,P.t=P.t,v2=v2,Yf=Yf)) #****
  } 
}
results<-kalman(para=para,Y=data,lik=lik,prev=prev,ahead=ahead) #**** 
results

# -2887.301 full sample

# Numerical Optimization  

# Lower bound & Upper bound 

low<-c(0.00001,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	 0.10,-0.15,-0.15,
	-0.15, 0.10,-0.15,
	-0.15,-0.15, 0.10,

	 -14.000,-14.000,-14.000,
	 
	 0,
	-0.99,0,
	-0.99,-0.99,0)
	
up<-c(0.99999,
	0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,
	 1.05,1.00,1.00,
	 1.00,1.00,1.00,
	 1.00,1.00,1.00,
	
	14.000,14.0000,14.0000,
	
	1.25,
	1.25,1.25,
	1.25,1.25,1.25)


#  nlminb

otim<-nlminb(para,kalman,Y=data,lik=lik,prev=prev,ahead=ahead,control = list(trace=1,iter.max=50000),lower=low,upper=up)

'
# Parameters after optimization

para<-c(0.07738452,
0.26726293,0.07443629,0.09040581,0.10480007,0.09931606,0.08645425,0.07852442,0.07207583,0.07269172,0.07909206,0.10295433,0.09261623,0.10050282,0.11161028,0.10670647,0.15070219,0.17327333,

0.99557203,0.02993130,-0.02087407,
-0.02531064,0.93776486,0.03662998,
0.03013432,0.02251098,0.83815771,

8.35291770,-1.44060409,-0.10671270,

0.30815585,
-0.04171767,0.61700758,
0.11418710,0.02219852,0.89692230)
'

lik<-TRUE
prev<-FALSE
results<-kalman(para=para,Y=data,lik=lik,prev=prev,ahead=ahead)
results

'
$objective
[1] -3184.239

$convergence
[1] 1

$iterations
[1] 147

$evaluations
function gradient 
     200     5507 

$message
[1] "function evaluation limit reached without convergence (9)"
'
