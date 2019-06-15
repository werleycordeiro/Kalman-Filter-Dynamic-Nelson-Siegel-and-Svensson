lyapunov<-function(N,phi,Q){
  matrix(solve(diag(N^2) - kronecker(phi,phi)) %*% matrix(Q,(N^2),1),N,N)
}
