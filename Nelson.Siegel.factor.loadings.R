Nelson.Siegel.factor.loadings <- function(l,m)
{
  column1 <- rep.int(1,length(m))
  column2 <- (1 - exp(-l * m))/(l * m)
  column3 <- column2 - exp(-l * m)
  
  lambmat <- cbind(column1,column2,column3)
  
  lambmat
}  
