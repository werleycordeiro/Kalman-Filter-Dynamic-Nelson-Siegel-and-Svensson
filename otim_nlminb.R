otim_nlminb = function(para,Y,lik,prev,ahead,matu,lower,upper){
  source("DNS-baseline.R")
  nlminb(para,kalman,Y=data,lik=lik,prev=prev,ahead=ahead,matu=matu,control = list(trace=1,iter.max=50000),lower=low,upper=up)
}
otim = otim_nlminb(para=para,Y=data,lik=lik,prev=prev,ahead=ahead,matu=matu,lower=low,upper=up)
return(otim)
