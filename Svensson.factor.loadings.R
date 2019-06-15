Svensson.factor.loadings <- function(para,m)
{
  l1   <- para[1]
  l2   <- para[2]  
  
  column1 <- rep.int(1,length(m))
  column2 <- (1 / (l1 * m)) * exp(-l1 * m) * (exp(l1 * m) - 1)
  column3 <- (1 / (l1 * m)) * exp(-l1 * m) * (-l1 * m + exp(l1 * m) - 1)
  column4 <- (1 / (l2 * m)) * exp(-l2 * m) * (-l2 * m + exp(l2 * m) - 1)
  
  lambmat <- cbind(column1,column2,column3,column4)
  
  lambmat
}  
