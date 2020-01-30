Estimate_hu <- function(rho, full_data, Z, mcmc0, mcmc1)
{
  library(mgcv)
  Npat <- nrow(full_data) #the number of patients
  Ncov <- ncol(Z) 
  id1 = which(Z[,1]==1)
  id0 = which(Z[,1]==0)
  N1 = length(id1)
  N0 = length(id0)
  cov = cbind(rep(1, Npat), Z[,-1])
  cov1 = Z[id1,]
  cov0 = cbind(rep(1, N0), Z[id0,-1])
  c1 = matrix(0, length(id1), Npat)
  c0 = matrix(0, length(id0), Npat)
  H=10
  temp0 = matrix(0, N0, N0)
  for (i in 1:Ncov)
  {
    temp0 = temp0 + outer(cov0[,i],cov0[,i],'-')^2
  }
  covariance0 = exp(-temp0) + 0.01*diag(N0)
  inv_covariance0 = solve(covariance0)
  temp1 = matrix(0, N1, N1)
  for (i in 1:Ncov)
  {
    temp1 = temp1 + outer(cov1[,i],cov1[,i],'-')^2
  }
  covariance1 = exp(-temp1) + 0.01*diag(N1)
  inv_covariance1 = solve(covariance1)
  for (k in 1:Npat)
  {
    covD = cov[k,]
    coeff1 = t(cov1)-covD
    c1[,k] = exp(-colSums(coeff1^2))
    coeff0 = t(cov0)-covD
    c0[,k] = exp(-colSums(coeff0^2))
  }
  var1 = 1.01-diag(t(c1)%*%inv_covariance1%*%c1)
  var0 = 1.01-diag(t(c0)%*%inv_covariance0%*%c0)
  partial_cov1 = t(c1)%*%inv_covariance1
  partial_cov0 = t(c0)%*%inv_covariance0
  
  u_range = seq(0,6,0.5)
  hu_result <- array(0, c(length(u_range), nsave, 2))
  for (u_index in 1:length(u_range))
    hu_result[u_index,,] = fbar_iter(mcmc0, mcmc1, u_range[u_index], rho, Npat,partial_cov0,partial_cov1, cov0, cov1, cov, var1, var0)
  return(hu_result)
}
