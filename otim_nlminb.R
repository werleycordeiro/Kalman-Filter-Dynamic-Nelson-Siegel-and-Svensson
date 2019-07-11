otim_nlminb = function(model,para,Y,lik,prev,ahead,matu,lower,upper){
  source(paste(model))
  nlminb(para,kalman,Y=data,lik=lik,prev=prev,ahead=ahead,matu=matu,control = list(trace=1,iter.max=50000),lower=low,upper=up)
}
otim = otim_nlminb(model=model,para=para,Y=data,lik=lik,prev=prev,ahead=ahead,matu=matu,lower=low,upper=up)
return(otim)
