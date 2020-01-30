Marginal_survival <- function(full_data, Z, mcmc0, mcmc1)
{
  Npat <- nrow(full_data) #the number of patients
  Ncov <- ncol(Z) 


  id1 = which(Z[,1]==1)
  id0 = which(Z[,1]==0)
  N1 = length(id1)
  N0 = length(id0)
  cov = cbind(rep(1, Npat), Z[,-1])
  cov1 = Z[id1,]
  cov0 = cbind(rep(1, N0), Z[id0,-1])
  
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
  
  
  tim = seq(0,10,0.3)
  fmean0 = array(NA, c(Npat,nsave, length(tim)))
  fmean1 = array(NA, c(Npat,nsave, length(tim)))
  for (k in 1:Npat)
  {
    #print(k)
    covD = cov[k,]
    fgrid0 = NULL
    fgrid1 = NULL
    for (i in 1:nsave)
    {
      fmean0[k,i,] <- 1-fmar_survival(tim, mcmc0$wh[,i],mcmc0$muh[,2,,i], mcmc0$betah2[,,i], mcmc0$Sigma[2,2,i], covD,cov0,inv_covariance0)
      fmean1[k,i,] <- 1-fmar_survival(tim, mcmc1$wh[,i],mcmc1$muh[,2,,i], mcmc1$betah2[,,i], mcmc1$Sigma[2,2,i], covD,cov1,inv_covariance1)
    }
  }
  tmp0 = matrix(0, nsave, length(tim))
  tmp1 = matrix(0, nsave, length(tim))
  for (i in 1:nsave) 
  {
    tmp0[i,]=apply(fmean0[,i,], 2, mean)
    tmp1[i,]=apply(fmean1[,i,], 2, mean)
  }
  fmean0_ave = apply(tmp0, 2, mean)
  fmean1_ave = apply(tmp1, 2, mean)
  fquantile0 = apply(tmp0, 2, function(x) quantile(x, c(0.025, 0.975)))
  fquantile1 = apply(tmp1, 2, function(x) quantile(x, c(0.025, 0.975)))
  return(list(fmean0_ave=fmean0_ave, fmean1_ave=fmean1_ave, fquantile0=fquantile0, fquantile1=fquantile1))
}