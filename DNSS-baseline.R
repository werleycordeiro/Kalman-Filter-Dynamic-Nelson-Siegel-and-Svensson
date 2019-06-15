# author: Werley Cordeiro
# werleycordeiro@gmail.com

# Data

data<-read.csv("https://www.dropbox.com/s/inpnlugzkddp42q/bonds.csv?dl=1",header = TRUE, sep = ";")
data<- as.matrix(data[,which(names(data)=="M3"):which(names(data)=="M120")])
datas<- seq(as.Date("1972/1/1"), by = "month", length.out = 348) # Parte 1: 1:116. Parte 2: 117:232. Parte 3: 233:348. 
require(xts)
data<- xts(data, order.by = datas)
data1<- xts(data, order.by = datas) #***

# Initial Parameters (See Dynamic-Nelson-Siegel to initial Parameters)  
# Lambda(2), H diag matrix (17x17), phi matrix (4x4), mu vector (4x1), Q matrix (4x4). Total: 49 parameters

para11<-c(0.077,0.077,
          0.26817668,0.07544778,0.09030891,0.10452337,0.09916127,0.08648250,0.07862852,0.07209767,0.07269271,0.07911056,0.10296477,0.09261803,0.10042012,0.11176658,0.10697570,0.15070754,0.17279775,
          
          0.99434649,0.02858860,-0.02207605,0.001,
          -0.02890870,0.93901662,0.03953932,0.001,
          0.02514169,0.02280013,0.841683100,0.001,
          0.00000001,0.00000001,0.000000001,0.986,
          
          8.35738231,-1.41976865,-0.15015063,1.78,
          
          0.30760631,
          -0.04519515,0.14206861,
          0.01694981,0.02563889,0.88259125,
          0.00000000,0.00000000,0.00000000,0.690)


prev<- FALSE # TRUE to Forecast
ahead<-12 # X-step ahead forecast
lik<-TRUE # # TRUE to return the value of the log-likelihood function. FALSE to return parameters.

# Kalman Filter function

kalman11<-function(para,Y,lik,prev,ahead){
  
m <- c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)  
M<-ahead#***

# Resize data if Forecast is on.
  if(prev){#***
    T <- nrow(Y)
    Yf<-Y
    Yf[(T-M+1):T,]<-NA
    Y<-Y[1:(T-M),]
    T <- nrow(Y)
  }else{
    T <- nrow(Y)
    Yf<-1
  }#***
  
pars<-list()
W <- ncol(Y)
N <- 4
  
# Create vectors and matrices
  pars$mu	<- matrix(NA,N,1)
  pars$phi	<- diag(N)
  pars$H	<- diag(ncol(Y))
  pars$Q	<- diag(N)
  
# Loading matrix
  
source("Svensson.factor.loadings.R")
pars$Z	<- Svensson.factor.loadings(para=para,m=m)
  
# Variance matrix of residuals
  for(i in 1:17){
    pars$H[i,i]<-para[i+2]
  }
  
H <- pars$H %*% pars$H
  
# Vector autoregressive coeffient matrix: VAR(1)
pars$phi[1,1] <- para[20]
pars$phi[1,2] <- para[21]
pars$phi[1,3] <- para[22]
pars$phi[1,4] <- para[23]
pars$phi[2,1] <- para[24]
pars$phi[2,2] <- para[25]
pars$phi[2,3] <- para[26]
pars$phi[2,4] <- para[27]
pars$phi[3,1] <- para[28]
pars$phi[3,2] <- para[29]
pars$phi[3,3] <- para[30]
pars$phi[3,4] <- para[31]
pars$phi[4,1] <- para[32]
pars$phi[4,2] <- para[33]
pars$phi[4,3] <- para[34]
pars$phi[4,4] <- para[35]
  
# Mean vector
pars$mu[1]<-para[36]
pars$mu[2]<-para[37]
pars$mu[3]<-para[38]
pars$mu[4]<-para[39]
  
# transition covariance matrix of residuals
pars$Q[1,1] <- para[40]
pars$Q[2,1] <- para[41]
pars$Q[2,2] <- para[42]
pars$Q[3,1] <- para[43]
pars$Q[3,2] <- para[44]
pars$Q[3,3] <- para[45]
pars$Q[4,1] <- para[46]
pars$Q[4,2] <- para[47]
pars$Q[4,3] <- para[48]
pars$Q[4,4] <- para[49]
  
Q <- pars$Q %*% t(pars$Q) 
  
v1   <- matrix(NA,T,W)			
v2   <- matrix(NA,T,W)
  
# Resize data if Forecast is on.
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
  
a.t[1, ]  <- pars$mu #pars$at0

# Start variance matrix
source("lyapunov.R")  
P.t[1, ,] <-lyapunov(N=N,phi=pars$phi,Q=Q) # Start variance matrix. pars$Pt0
  
# Initial loglikelihood	  
logLik <- - 0.5 * T * ncol(Y) * log(2 * pi)

# Kalman Filter and log-likelihood
source("Kfilter.R")  
Kfilter(logLik=logLik,N=N,T=T,Y=Y,Z=pars$Z,a.t=a.t,P.t=P.t,H=H,a.tt=a.tt,P.tt=P.tt,v2=v2,phi=pars$phi,mu=pars$mu,Q=Q,prev=prev,M=M,Yf=Yf,lik=lik)
}

results11<-kalman11(para=para11,Y=data,lik=lik,prev=prev,ahead=ahead) 
results11

# -1839.368


#  Otimização


#---------------#
#  Nelder-Mead  #
#---------------#
'
otim <- optim(par=para,fn=kalman1,Y=data,lik=lik,gr=NULL,method="Nelder-Mead",
control=list(trace=2,REPORT=2,maxit=10000000),hessian = FALSE)
'

#---------#
# LB & UB #
#---------#

low11<-c(0.0000001,0.0000001,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
         0.10,-0.15,-0.15,-0.15,
         -0.15, 0.10,-0.15,-0.15,
         -0.15,-0.15, 0.10,-0.15,
         -0.15,-0.15,-0.15, 0.10,
         
         1.000,-5.000,-2.000,-5.000,
         
         0.001,
         -1.25,0.001,
         -1.25,-1.25,0.001,
         -1.25,-1.25,-1.25,0.001)

up11<-c(0.99999,0.99999,
        0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,
        1.05,0.50,0.50,0.50,
        0.50,1.05,0.50,0.50,
        0.50,0.50,1.05,0.50,
        0.50,0.50,0.50,1.05,
        
        10,14,14,10,
        
        1.25,
        1.00,1.25,
        1.00,1.00,1.25, 
        1.00,1.00,1.00,1.25)

#----------#
#  nlminb  #
#----------#
'
otim1<-nlminb(para11,kalman11,Y=data,lik=lik,gradient = NULL, hessian = NULL,control = list(trace=1,iter.max=50000,eval.max=50000,rel.tol=1e-5),
lower=low11,upper=up11)
'
#-----------------------#
#  optim/optimParallel  #
#-----------------------#
'
checkin<-Sys.time()

require(optimParallel)
cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl) #***************************

otim <- optim(par=para1,fn=kalman1,Y=data,lik=lik,gr=NULL,method="L-BFGS-B", #********************* "optimParallel"
lower=low,upper=up,control=list(trace=3,REPORT=2,maxit=2000),hessian = FALSE)

checkout<-Sys.time()
'
para11<-c(0.104467920,0.049441788,
          0.223430801,0.000000000,0.094984312,0.103427243,0.077322872,0.064905752,0.068974548,0.08029856,0.070889403,0.064934593,0.085058637,0.065476346,0.100203275,0.114681066,0.096374289,0.108923682,0.172100969,
          
          1.001605218,0.032684878,-0.027546959,-0.009417855,
          -0.038057180,0.927345689,0.080668263,0.027995407,
          0.036720125,0.064640757,0.789410433,-0.028201668,
          -0.010575874,0.023509296,0.037788015,0.913207235,
          
          7.946739213,-0.663136893,0.405595470,-0.300647232,
          
          0.338044756,
          -0.152405470,0.704075030,
          0.425986914,-0.369501881,1.152326423,
          -0.302120250,0.344675961,-0.373388285,0.693240448)

lik<-TRUE
prev<-FALSE
results11<-kalman11(para=para11,Y=data,lik=lik,prev=prev,ahead) 
results11

# -3678.985 full sample
