#-----------------#
# 	DNS-baseline  #
#-----------------#
list.of.packages <- c("optimParallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(optimParallel)
#------------------#
# 		Dados      #
#------------------#

#setwd("C:\\Users\\werle\\Dropbox\\Mestrado_UFSC\\4_Semestre\\DNS\\DNS-baseline")
#data <- read.csv("bonds.csv",header = TRUE, sep = ";")
data<-read.csv("https://www.dropbox.com/s/inpnlugzkddp42q/bonds.csv?dl=1",header = TRUE, sep = ";")#***
require(xts)
data<- as.matrix(data[,which(names(data)=="M3"):which(names(data)=="M120")])
datas<- seq(as.Date("1972/1/1"), by = "month", length.out = 348)
data<- xts(data, order.by = datas)
data1<- xts(data, order.by = datas) #***


#------------------------------#
#    Parâmetros Iniciais	   #
#------------------------------#

# 36
# Lambda, Matriz H (17x17), Matriz Phi(3x3),Vetor mu (4x1), Matriz Q (3x3), 
 
para<-c(0.0609,

0.14170940,0.07289485,0.11492339,0.11120008,0.09055795,0.07672075,0.07222108,0.07076431,0.07012891,0.07267366,0.10624206,0.09029621,0.10374527,0.09801215,0.09122014,0.11794190,0.13354418,

 0.99010443,0.02496842,-0.002294319,
-0.02812401,0.94256154, 0.028699387,
 0.05178493,0.01247332, 0.788078795, 

8.345444,-1.572442,0.2029919,  

 0.3408764,
-0.07882772,0.62661018,
-0.21351036,-0.00425989,1.08802059)

#------------------#
# Filtro de Kalman #
#------------------#

prev<- FALSE #***
ahead<-12 #***
lik <- TRUE # Porque quero analisar somente o valor do loglikelihood.

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
  
  M<-ahead#***
  
	  if(prev){#***
	  T <- nrow(Y)
	  Yf<-Y
	  Yf[(T-M+1):T,]<-NA
	  Y<-Y[1:(T-M),]
	  T <- nrow(Y)
	  }else{
	  T <- nrow(Y)
	  }#***
	  
  pars<-list()
  W <- ncol(Y)
  N <- 3
  
  pars$mu	<- matrix(NA,N,1)
  pars$phi	<- diag(N)
  pars$H	<- diag(ncol(Y))
  pars$Q	<- diag(N)
  pars$Z	<- Nelson.Siegel.factor.loadings(l)
  
  # disturbance covariance matrix 
  for(i in 1:17){
  pars$H[i,i]<-para[1 + i]
  }

  H <- pars$H^2
  
  # Matriz de coeficientes do VAR.
  pars$phi[1,1] <- para[19]
  pars$phi[1,2] <- para[20]
  pars$phi[1,3] <- para[21]
  pars$phi[2,1] <- para[22]
  pars$phi[2,2] <- para[23]
  pars$phi[2,3] <- para[24]
  pars$phi[3,1] <- para[25]
  pars$phi[3,2] <- para[26]
  pars$phi[3,3] <- para[27]
  
  # 3 mean state vector
  pars$mu[1]<-para[28]
  pars$mu[2]<-para[29]
  pars$mu[3]<-para[30]
  
  # transition covariance matrix of residuals
  pars$Q[1,1] <- para[31]
  pars$Q[2,1] <- para[32]
  pars$Q[2,2] <- para[33]
  pars$Q[3,1] <- para[34]
  pars$Q[3,2] <- para[35]
  pars$Q[3,3] <- para[36]
  
  Q <- pars$Q %*% t(pars$Q) 
  
  v2   <- matrix(NA,T,W)			
  
  if(prev){#***
	  a.tt <- matrix(NA, (T+M), N)
	  a.t  <- matrix(NA, (T+M+1), N) # caso prev=TRUE, sempre será dim(a.t)[1]=348
	  P.tt <- array(NA, c((T+M), N, N))
	  P.t  <- array(NA, c((T+M+1), N, N))
  }else{
	  a.tt <- matrix(NA, T, N)
	  a.t  <- matrix(NA, (T+1), N)
	  P.tt <- array(NA, c(T, N, N))
	  P.t  <- array(NA, c((T+1), N, N))
  }#***
  
  a.t[1, ]  <- pars$mu #pars$at0
  P.t[1, ,] <- matrix(solve(diag(N^2) - kronecker(pars$phi,pars$phi)) %*% matrix(Q,(N^2),1),N,N)  # pars$Pt0 #
  
  logLik <- - 0.5 * T * ncol(Y) * log(2 * pi)
  
  # Filtro
  for (t in 1:T) 
  {
    # Log-likelihood
	v <- (as.numeric(Y[t, ])) - pars$Z %*% a.t[t, ] # erro de previsão modificado
    F <- pars$Z %*% P.t[t, ,] %*% t(pars$Z) + H # matriz var-cov do erro de previsão
		if(det(F)<=1e-30 || is.na(det(F)) || is.nan(det(F)) || is.infinite(det(F))){
			logLik<- -1000000000000000; break
		}else{
	F.inv  <- solve(F)
    logLik <- logLik - 0.5 * (log(det(F)) + t(v) %*% F.inv %*% v)
	# Update
    a.tt[t, ]   <- a.t[t, ] +  P.t[t, , ] %*% t(pars$Z) %*% F.inv %*% v
	P.tt[t, , ] <- P.t[t, , ] - P.t[t, , ] %*% t(pars$Z) %*% F.inv %*% pars$Z %*% P.t[t, , ]
    v2[t, ] <- (as.numeric(Y[t, ])) - pars$Z %*% a.tt[t, ]
	# Predict
	a.t[t + 1, ]  <- pars$phi %*% a.tt[t, ] + (diag(N) - pars$phi) %*% pars$mu  
	P.t[t + 1, ,] <- pars$phi %*% P.tt[t, ,] %*% t(pars$phi) + Q
	}
		if(prev) # Forecast
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

#--------------#
#  Otimização  #
#--------------#

#-----------#
#  LB & UB  #
#-----------#

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

#----------#
#  nlminb  # full sample
#----------#
' 
otim<-nlminb(para,kalman,Y=data,lik=lik,prev=prev,ahead=ahead,control = list(trace=1,iter.max=50000),lower=low,upper=up)
'
para<-c(0.07738452,
0.26726293,0.07443629,0.09040581,0.10480007,0.09931606,0.08645425,0.07852442,0.07207583,0.07269172,0.07909206,0.10295433,0.09261623,0.10050282,0.11161028,0.10670647,0.15070219,0.17327333,

0.99557203,0.02993130,-0.02087407,
-0.02531064,0.93776486,0.03662998,
0.03013432,0.02251098,0.83815771,

8.35291770,-1.44060409,-0.10671270,

0.30815585,
-0.04171767,0.61700758,
0.11418710,0.02219852,0.89692230)

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

#---------------------------------#
# Gráfico p/ teste na full sample #
#---------------------------------#

#ts.plot(as.numeric(data1[337,]),type="p",ylim=c(4,7))
#lines(as.numeric(results$Yf[337,]))

#----------#
# Previsão #
#----------#

# Dividi a amostra em quatro partições: P1:[1:261,] até [21:281,]  #21 janelas
#										P2:[22:282,] até [42:302,] #21 janelas
#										P3:[43:303,] até [63:323,] #21 janelas
#										P4:[64:324,] até [88:348,] #25 janelas. Total = 21 + 21 + 21 + 25 = 88 janelas
# Em cada janela de cada partição, isto é, na partição 1, janela [1:261,], [2:262,], [3:263,], ..., os parâmetros foram estimados
# baseados nos dados até t-12. Por exemplo: em P1 [1:249,], os parâmetros são estimados e então as previsões 12 passos a frente são feitas.  

#-------------#
# Partição 1  #
#-------------#
checkin<-Sys.time()
mse<-array(NA,c(12,17,21)) # mudar último elemento dependendo do número de previsões com otimização
colnames(mse)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120")
rownames(mse)<-c("1 month","2 months","3 months","4 months","5 months","6 months",
                  "7 months","8 months","9 months","10 months","11 months","12 months")
				  
for(j in 1:21){
data2<-data[j:(260+j),]
lik<-TRUE
prev<- FALSE
otim<-nlminb(para,kalman,Y=data2,lik=lik,prev=prev,ahead=ahead,control = list(trace=1,iter.max=500),lower=low,upper=up)
lik<-FALSE
prev<-TRUE
ahead<-12
results<-kalman(para=otim$par,Y=data2,lik=lik,prev=prev,ahead=ahead)
	for(i in 1:12)
	{
	  mse[i,,j]<-(results$Yf[i+(261-12),]-data2[i+(261-12),])^2 # falta calcular a média e tirar a raiz quadrada. 13 pq ahead+1
	}
print(c("Finalizado iteração & previsão DNS-base---->",j))
}
mse1.DNS.base<-mse
save(mse1.DNS.base,file="mse1.DNS.base.rda")
checkout<-Sys.time()
checkout-checkin
#-------------#
# Partição 2  #
#-------------#
checkin<-Sys.time()
mse<-array(NA,c(12,17,21)) # mudar último elemento dependendo do número de previsões com otimização
colnames(mse)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120")
rownames(mse)<-c("1 month","2 months","3 months","4 months","5 months","6 months",
                  "7 months","8 months","9 months","10 months","11 months","12 months")
				  
for(j in 1:21){
data2<-data[(21+j):(281+j),]
lik<-TRUE
prev<- FALSE
otim<-nlminb(para,kalman,Y=data2,lik=lik,prev=prev,ahead=ahead,control = list(trace=1,iter.max=500),lower=low,upper=up)
lik<-FALSE
prev<-TRUE
ahead<-12
results<-kalman(para=otim$par,Y=data2,lik=lik,prev=prev,ahead=ahead)
	for(i in 1:12)
	{
	  mse[i,,j]<-(results$Yf[i+(261-12),]-data2[i+(261-12),])^2 # falta calcular a média e tirar a raiz quadrada. 13 pq ahead+1
	}
print(c("Finalizado iteração & previsão DNS-base---->",21+j))
}
mse2.DNS.base<-mse
save(mse2.DNS.base,file="mse2.DNS.base.rda")
checkout<-Sys.time()
checkout-checkin
#-------------#
# Partição 3  #
#-------------#
checkin<-Sys.time()
mse<-array(NA,c(12,17,21)) # mudar último elemento dependendo do número de previsões com otimização
colnames(mse)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120")
rownames(mse)<-c("1 month","2 months","3 months","4 months","5 months","6 months",
                  "7 months","8 months","9 months","10 months","11 months","12 months")
				  
for(j in 1:21){
data2<-data[(42+j):(302+j),]
lik<-TRUE
prev<- FALSE
otim<-nlminb(para,kalman,Y=data2,lik=lik,prev=prev,ahead=ahead,control = list(trace=1,iter.max=500),lower=low,upper=up)
lik<-FALSE
prev<-TRUE
ahead<-12
results<-kalman(para=otim$par,Y=data2,lik=lik,prev=prev,ahead=ahead)
	for(i in 1:12)
	{
	  mse[i,,j]<-(results$Yf[i+(261-12),]-data2[i+(261-12),])^2 # falta calcular a média e tirar a raiz quadrada. 13 pq ahead+1
	}
print(c("Finalizado iteração & previsão DNS-base---->",42+j))
}
mse3.DNS.base<-mse
save(mse3.DNS.base,file="mse3.DNS.base.rda")
checkout<-Sys.time()
checkout-checkin
#-------------#
# Partição 4  #
#-------------#
checkin<-Sys.time()
mse<-array(NA,c(12,17,25)) # mudar último elemento dependendo do número de previsões com otimização
colnames(mse)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120")
rownames(mse)<-c("1 month","2 months","3 months","4 months","5 months","6 months",
                  "7 months","8 months","9 months","10 months","11 months","12 months")
				  
for(j in 1:25){
data2<-data[(63+j):(323+j),]
lik<-TRUE
prev<- FALSE
otim<-nlminb(para,kalman,Y=data2,lik=lik,prev=prev,ahead=ahead,control = list(trace=1,iter.max=500),lower=low,upper=up)
lik<-FALSE
prev<-TRUE
ahead<-12
results<-kalman(para=otim$par,Y=data2,lik=lik,prev=prev,ahead=ahead)
	for(i in 1:12)
	{
	  mse[i,,j]<-(results$Yf[i+(261-12),]-data2[i+(261-12),])^2 # falta somar, dividir por R e tirar a raiz quadrada. 13 pq ahead+1
	}
print(c("Finalizado iteração & previsão DNS-base---->",63+j))
}
mse4.DNS.base<-mse
save(mse4.DNS.base,file="mse4.DNS.base.rda")
checkout<-Sys.time()
checkout-checkin

#-----------------------------------#
# Analisar resultados das previsões #
#-----------------------------------#
'
setwd("C:\\Users\\werle\\Dropbox\\Mestrado_UFSC\\4_Semestre\\DNS\\DNS-baseline\\Previsão")

load(file = "mse1.DNS.base.rda")
load(file = "mse2.DNS.base.rda")
load(file = "mse3.DNS.base.rda")
load(file = "mse4.DNS.base.rda")
library(abind)
mse<-abind(mse1.DNS.base, mse2.DNS.base,mse3.DNS.base,mse4.DNS.base, along = 3)
dim(mse)

rmspe.dns<-matrix(NA,12,17)
colnames(rmspe.dns)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120")
rownames(rmspe.dns)<-c("1 month","2 months","3 months","4 months","5 months","6 months",
                  "7 months","8 months","9 months","10 months","11 months","12 months")

for(j in 1:12){
	for(i in 1:17){
		rmspe.dns[j,i]<-(sqrt(mean(mse[j,i,])))*100
	}
}
rmspe.dns
'
#----#
# A1 #
#----#
setwd("C:\\Users\\werle\\Dropbox\\Mestrado_UFSC\\4_Semestre\\DNS\\DNS-baseline\\Previsão")
load(file = "mse1.DNS.base.rda")
mse<-mse1.DNS.base

rmspe.dns<-matrix(NA,12,17)
colnames(rmspe.dns)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120")
rownames(rmspe.dns)<-c("1 month","2 months","3 months","4 months","5 months","6 months",
                  "7 months","8 months","9 months","10 months","11 months","12 months")
for(j in 1:12){
	for(i in 1:17){
		rmspe.dns[j,i]<-(sqrt(mean(mse[j,i,])))*100
	}
}
rmspe.dns

#----#
# A2 #
#----#
setwd("C:\\Users\\werle\\Dropbox\\Mestrado_UFSC\\4_Semestre\\DNS\\DNS-baseline\\Previsão")
load(file = "mse2.DNS.base.rda")
mse<-mse2.DNS.base

rmspe.dns<-matrix(NA,12,17)
colnames(rmspe.dns)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120")
rownames(rmspe.dns)<-c("1 month","2 months","3 months","4 months","5 months","6 months",
                  "7 months","8 months","9 months","10 months","11 months","12 months")
for(j in 1:12){
	for(i in 1:17){
		rmspe.dns[j,i]<-(sqrt(mean(mse[j,i,])))*100
	}
}
rmspe.dns

#----#
# A3 #
#----#
setwd("C:\\Users\\werle\\Dropbox\\Mestrado_UFSC\\4_Semestre\\DNS\\DNS-baseline\\Previsão")
load(file = "mse3.DNS.base.rda")
mse<-mse3.DNS.base

rmspe.dns<-matrix(NA,12,17)
colnames(rmspe.dns)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120")
rownames(rmspe.dns)<-c("1 month","2 months","3 months","4 months","5 months","6 months",
                  "7 months","8 months","9 months","10 months","11 months","12 months")
for(j in 1:12){
	for(i in 1:17){
		rmspe.dns[j,i]<-(sqrt(mean(mse[j,i,])))*100
	}
}
rmspe.dns

#----#
# A4 #
#----#
setwd("C:\\Users\\werle\\Dropbox\\Mestrado_UFSC\\4_Semestre\\DNS\\DNS-baseline\\Previsão")
load(file = "mse4.DNS.base.rda")
mse<-mse4.DNS.base

rmspe.dns<-matrix(NA,12,17)
colnames(rmspe.dns)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120")
rownames(rmspe.dns)<-c("1 month","2 months","3 months","4 months","5 months","6 months",
                  "7 months","8 months","9 months","10 months","11 months","12 months")
for(j in 1:12){
	for(i in 1:17){
		rmspe.dns[j,i]<-(sqrt(mean(mse[j,i,])))*100
	}
}
rmspe.dns

#------------------#
# Tabela Previsões #
#------------------#
'
m1<-matrix(NA,10,17)
colnames(m1)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120")
rownames(m1)<-c("DNS-base","DNS-TVL","DNS-GARCH","DNS-GARCH-TVL","DNSS","DNSS-TVL1","DNSS-TVL2","DNSS-GARCH","DNSS-GARCH-TVL1","DNSS-GARCH-TVL2")
m3<-m1
m6<-m1
m12<-m1

m1[1,]<-rmspe.dns[1,]
m3[1,]<-rmspe.dns[3,]
m6[1,]<-rmspe.dns[6,]
m12[1,]<-rmspe.dns[12,]

setwd("C:\\Users\\werle\\Dropbox\\Mestrado_UFSC\\4_Semestre\\Tabela Previsões")
save(m1,file="m1.rda")
save(m3,file="m3.rda")
save(m6,file="m6.rda")
save(m12,file="m12.rda")
'
#----#
# A1 #
#----#
A1m1<-matrix(NA,10,17)
colnames(A1m1)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120")
rownames(A1m1)<-c("DNS-base","DNS-TVL","DNS-GARCH","DNS-GARCH-TVL","DNSS","DNSS-TVL1","DNSS-TVL2","DNSS-GARCH","DNSS-GARCH-TVL1","DNSS-GARCH-TVL2")
A1m3<-A1m1
A1m6<-A1m1
A1m12<-A1m1

A1m1[1,]<-rmspe.dns[1,]
A1m3[1,]<-rmspe.dns[3,]
A1m6[1,]<-rmspe.dns[6,]
A1m12[1,]<-rmspe.dns[12,]

setwd("C:\\Users\\werle\\Dropbox\\Mestrado_UFSC\\4_Semestre\\Tabela Previsões")
save(A1m1,file="A1m1.rda")
save(A1m3,file="A1m3.rda")
save(A1m6,file="A1m6.rda")
save(A1m12,file="A1m12.rda")

#----#
# A2 #
#----#
A2m1<-matrix(NA,10,17)
colnames(A2m1)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120")
rownames(A2m1)<-c("DNS-base","DNS-TVL","DNS-GARCH","DNS-GARCH-TVL","DNSS","DNSS-TVL1","DNSS-TVL2","DNSS-GARCH","DNSS-GARCH-TVL1","DNSS-GARCH-TVL2")
A2m3<-A2m1
A2m6<-A2m1
A2m12<-A2m1

A2m1[1,]<-rmspe.dns[1,]
A2m3[1,]<-rmspe.dns[3,]
A2m6[1,]<-rmspe.dns[6,]
A2m12[1,]<-rmspe.dns[12,]

setwd("C:\\Users\\werle\\Dropbox\\Mestrado_UFSC\\4_Semestre\\Tabela Previsões")
save(A2m1,file="A2m1.rda")
save(A2m3,file="A2m3.rda")
save(A2m6,file="A2m6.rda")
save(A2m12,file="A2m12.rda")

#----#
# A3 #
#----#
A3m1<-matrix(NA,10,17)
colnames(A3m1)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120")
rownames(A3m1)<-c("DNS-base","DNS-TVL","DNS-GARCH","DNS-GARCH-TVL","DNSS","DNSS-TVL1","DNSS-TVL2","DNSS-GARCH","DNSS-GARCH-TVL1","DNSS-GARCH-TVL2")
A3m3<-A3m1
A3m6<-A3m1
A3m12<-A3m1

A3m1[1,]<-rmspe.dns[1,]
A3m3[1,]<-rmspe.dns[3,]
A3m6[1,]<-rmspe.dns[6,]
A3m12[1,]<-rmspe.dns[12,]

setwd("C:\\Users\\werle\\Dropbox\\Mestrado_UFSC\\4_Semestre\\Tabela Previsões")
save(A3m1,file="A3m1.rda")
save(A3m3,file="A3m3.rda")
save(A3m6,file="A3m6.rda")
save(A3m12,file="A3m12.rda")

#----#
# A4 #
#----#
A4m1<-matrix(NA,10,17)
colnames(A4m1)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120")
rownames(A4m1)<-c("DNS-base","DNS-TVL","DNS-GARCH","DNS-GARCH-TVL","DNSS","DNSS-TVL1","DNSS-TVL2","DNSS-GARCH","DNSS-GARCH-TVL1","DNSS-GARCH-TVL2")
A4m3<-A4m1
A4m6<-A4m1
A4m12<-A4m1

A4m1[1,]<-rmspe.dns[1,]
A4m3[1,]<-rmspe.dns[3,]
A4m6[1,]<-rmspe.dns[6,]
A4m12[1,]<-rmspe.dns[12,]

setwd("C:\\Users\\werle\\Dropbox\\Mestrado_UFSC\\4_Semestre\\Tabela Previsões")
save(A4m1,file="A4m1.rda")
save(A4m3,file="A4m3.rda")
save(A4m6,file="A4m6.rda")
save(A4m12,file="A4m12.rda")

#----------------------#
# optim/optimParallel  #
#----------------------#
'	
require(optimParallel)
cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl) 
checkin<-Sys.time()
otim <- optim(par=para,fn=kalman,Y=data,lik=lik,prev=prev,ahead=ahead,gr=NULL,method="L-BFGS-B",
                lower=low,upper=up,control=list(trace=3,REPORT=2,maxit=50000),hessian = FALSE)
checkout<-Sys.time()
'

#--------#
#  StDev #
#--------#
'
hessian<-optimHess(par=otim$par, fn=kalman,Y=data,lik=lik)
lik<-TRUE
sqrt(diag(solve(hessian)))


sd<-c(0.002102698,
0.014004974,0.008946594,0.004508075,0.004639367,0.004414822,0.004051509,0.003824602,0.003462985,0.003594687,0.003777324,0.004377243,0.004227602,0.004540631,0.005085867,0.005399538,0.006857300,0.007624306,

0.008758588,0.008999350,0.010830216,
0.014377874,0.017882707,0.021397333,
0.028195251,0.026758894,0.031920421,

1.484348252,0.475985220,0.453283952,

0.013722009,
0.036544496,0.024540563,
0.061871426,0.056818688,0.048222808)
'
#----------#
#  lbfgsb3 #
#----------#
'
require(lbfgsb3)
otim1 <- lbfgsb3(prm=para,fn=kalman,Y=data,lik=lik,prev=prev,ahead=ahead,gr=NULL,lower=low,upper=up,control=list(trace=3,maxit=2000))
'
#-----------------#
# Erros filtrados #
#-----------------#

# Média Obs.: Calcular sobre os parâmetros ótimos
results<-kalman(para=para,Y=data,lik=FALSE)
md<-matrix(NA,20,2)
rownames(md)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120","Média","Mediana","Menor")
colnames(md)<-c("Média","Dp")
for(i in 1:17){md[i,1]<-mean(results$v2[,i]);
				md[i,2]<-sd(results$v2[,i]) ;
				md[18,1]<-mean(md[(1:17),1]);
				md[18,2]<-mean(md[(1:17),2]);
				md[19,1]<-median(md[(1:17),1]);
				md[19,2]<-median(md[(1:17),2])
				}
md<-md*100
md
'
              Média        Dp
M3      -12.4913557 22.303389
M6       -1.2482456  4.781173
M9        0.5398494  8.148333
M12       1.3238106  9.934143
M15       3.7108978  8.759942
M18       3.5741951  7.232261
M21       3.2078520  6.429982
M24      -1.4269662  6.327344
M30      -2.6772523  5.983386
M36      -3.2680055  6.624304
M48      -1.8690450  9.671982
M60      -3.2962585  7.968476
M72       1.9686366  9.015500
M84       0.6921133 10.152027
M96       3.4884691  9.111029
M108      4.1969236 13.496832
M120     -1.3030553 16.334162
Média    -0.2869080  9.545545
Mediana   0.5398494  8.759942
Menor            NA        NA

'
# AIC
2*36-2*3184.239


#saveRDS(para.vet,"para.vet.rds")
#para.vet<-readRDS("para.vet.rds") # Usar no .Rmd


#saveRDS(datacomp,"datacomp2.rds")
#datacomp<-readRDS("datacomp.rds") # Usar no .Rmd

#saveRDS(data,"data2.rds")
#data<-readRDS("data.rds") # Usar no .Rmd


# Gráfico 1
par(oma = c(4, 1, 1, 1))
par(mfrow=c(3,4))

ts.plot(datacomp[325,],ylim=range(4:5),main="Previsão 1 : Jan-1999",xlab="Maturidades",ylab="Taxas")
lines(data[325,],col="blue",lty=2)
ts.plot(datacomp[326,],ylim=range(4:6),main="Previsão 2 : Fev-1999",xlab="Maturidades",ylab="Taxas")
lines(data[326,],col="blue",lty=2)
ts.plot(datacomp[327,],ylim=range(4:6),main="Previsão 3 : Mar-1999",xlab="Maturidades",ylab="Taxas")
lines(data[327,],col="blue",lty=2)
ts.plot(datacomp[328,],ylim=range(4:6),main="Previsão 4 : Abr-1999",xlab="Maturidades",ylab="Taxas")
lines(data[328,],col="blue",lty=2)
ts.plot(datacomp[329,],ylim=range(4:6),main="Previsão 5 : Mai-1999",xlab="Maturidades",ylab="Taxas")
lines(data[329,],col="blue",lty=2)
ts.plot(datacomp[330,],ylim=range(4:6),main="Previsão 6 : Jun-1999",xlab="Maturidades",ylab="Taxas")
lines(data[330,],col="blue",lty=2)

ts.plot(datacomp[331,],ylim=range(4:7),main="Previsão 7 : Jul-1999",xlab="Maturidades",ylab="Taxas")
lines(data[331,],col="blue",lty=2)
ts.plot(datacomp[332,],ylim=range(4:7),main="Previsão 8 : Ago-1999",xlab="Maturidades",ylab="Taxas")
lines(data[332,],col="blue",lty=2)
ts.plot(datacomp[333,],ylim=range(4:7),main="Previsão 9 : Set-1999",xlab="Maturidades",ylab="Taxas")
lines(data[333,],col="blue",lty=2)
ts.plot(datacomp[334,],ylim=range(4:7),main="Previsão 10 : Out-1999",xlab="Maturidades",ylab="Taxas")
lines(data[334,],col="blue",lty=2)
ts.plot(datacomp[335,],ylim=range(4:7),main="Previsão 11 : Nov-1999",xlab="Maturidades",ylab="Taxas")
lines(data[335,],col="blue",lty=2)
ts.plot(datacomp[336,],ylim=range(4:7),main="Previsão 12 : Dez-1999",xlab="Maturidades",ylab="Taxas")
lines(data[336,],col="blue",lty=2)


title(" ",line=-1, outer = TRUE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Observado", "Previsto"), xpd = TRUE, horiz = TRUE, inset = c(0, 
                                                                                 0), bty = "n", col = c("black", "blue"),lty=1:2, cex = 1)

# Root Mean Squared Error

rmse<-matrix(NA,12,17)
colnames(rmse)<-c("M3","M6","M9","M12","M15","M18","M21","M24","M30","M36","M48","M60","M72","M84","M96","M108","M120")
rownames(rmse)<-c("1 month","2 months","3 months","4 months","5 months","6 months",
                  "7 months","8 months","9 months","10 months","11 months","12 months")
for(i in 325:336)
{
  rmse[(i-324),]<-(sqrt((datacomp[i,]-data[i,])^2))
}
rmse

#saveRDS(rmse,"rmse1.rds") 1 - SEM CRISE
rmse<-readRDS("rmse1.rds") # Usar no .Rmd 

DNS<-matrix(NA,17,1)
for(i in 2:18){
DNS[(i-1),1]<-exp(0.08/sqrt(para[i]))/(1+exp(0.08/sqrt(para[i])))
}
DNS
