main_mcmc <- function(Niter, burn.in, lag, full_data, Z)
{
  Npat <- nrow(full_data) #the number of patients
  Ncov <- ncol(Z) 
  #global variable
  H = 10  #the upper limit for stick breaking of DP. 
  lambda0 = 4
  lambda1 = 1 #prior for M
  lambda2 = 1
  
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
  
  var1_betah = solve(t(cov1)%*%inv_covariance1%*%cov1+diag(Ncov))
  var0_betah = solve(t(cov0)%*%inv_covariance0%*%cov0+diag(Ncov))
  
  partial_cov0 = t(cov0)%*%inv_covariance0
  partial_cov1 = t(cov1)%*%inv_covariance1
  
  #prior for \Sigma
  Phi1 = cov(full_data[id1,1:2])
  Phi0 = cov(full_data[id0,1:2])
  
  
  mcmc1 <- NULL
  mcmc1$M <- rep(NA, Niter)
  mcmc1$Sigma <- array(NA, c(2,2,Niter))
  mcmc1$muh <- array(NA, c(H, 2, N1, Niter))
  mcmc1$betah1 <- array(NA, c(H, Ncov, Niter))
  mcmc1$betah2 <- array(NA, c(H, Ncov, Niter))
  mcmc1$r <- array(NA, c(N1, Niter))
  mcmc1$wh <- array(NA, c(H, Niter))
  mcmc1$imputey <- array(NA, c(N1,2, Niter))
  
  set.seed(1)
  ##Initialize 
  initial <- init(full_data[id1,1:2], N1, H)
  mcmc1$M[1] = initial$M
  mcmc1$Sigma[,,1] = initial$Sigma
  for (i in 1:N1) mcmc1$muh[,,i,1] = initial$mh
  mcmc1$wh[,1] = initial$wh
  mcmc1$r[,1] = initial$r
  lfit1 = survreg(Surv(full_data[id1,1],full_data[id1,3])~cov1,dist="gaussian")
  lfit2 = survreg(Surv(full_data[id1,2],full_data[id1,4])~cov1,dist="gaussian")
  mcmc1$betah1[,,1] = matrix(rep(lfit1$coefficients[-2], H),H,Ncov, byrow=T)
  mcmc1$betah2[,,1] = matrix(rep(lfit2$coefficients[-2], H),H,Ncov, byrow=T)
  mcmc1$imputey[,,1] = full_data[id1,1:2]
  prior_betah1 = mcmc1$betah1[1,,1]
  prior_betah2 = mcmc1$betah2[1,,1]
  
  print("MCMC Iteration under treatment:")
  for (iter in 2:Niter)
  {
    print(iter)
    tmp = update_r(mcmc1$Sigma[,,iter-1], mcmc1$muh[,,,iter-1],mcmc1$wh[,iter-1], full_data[id1,], mcmc1$r[,iter-1], N1)
    mcmc1$r[,iter] = tmp$r
    mcmc1$imputey[,,iter] = tmp$impute_y
    tmp2 = update_wh_and_M(mcmc1$r[,iter], mcmc1$M[iter-1],lambda1, lambda2)
    mcmc1$M[iter] = tmp2$M
    mcmc1$wh[,iter] = tmp2$wh
    mcmc1$muh[,,,iter] = update_muh(mcmc1$wh[,iter-1],mcmc1$r[,iter], mcmc1$Sigma[,,iter-1], mcmc1$muh[,,,iter-1], mcmc1$imputey[,,iter], mcmc1$betah1[,,iter-1],mcmc1$betah2[,,iter-1], N1, cov1,covariance1, inv_covariance1,partial_cov1)
    mcmc1$betah1[,,iter] = update_betah1(mcmc1$muh[,,,iter], mcmc1$r[,iter],prior_betah1,var1_betah,partial_cov1, Ncov)
    mcmc1$betah2[,,iter] = update_betah2(mcmc1$muh[,,,iter], mcmc1$r[,iter],prior_betah2,var1_betah,partial_cov1,Ncov)
    mcmc1$Sigma[,,iter] = update_Sigma(mcmc1$r[,iter], mcmc1$imputey[,,iter], mcmc1$muh[,,,iter], N1, Phi1,lambda0)
  }
  
  mcmc0 <- NULL
  mcmc0$M <- rep(NA, Niter)
  mcmc0$Sigma <- array(NA, c(2,2,Niter))
  mcmc0$muh <- array(NA, c(H, 2, N0, Niter))
  mcmc0$betah1 <- array(NA, c(H, Ncov, Niter))
  mcmc0$betah2 <- array(NA, c(H, Ncov, Niter))
  mcmc0$r <- array(NA, c(N0, Niter))
  mcmc0$wh <- array(NA, c(H, Niter))
  mcmc0$imputey <- array(NA, c(N0,2, Niter))
  
  set.seed(1)
  ##Initialize 
  initial <- init(full_data[id0,1:2], N0, H)
  mcmc0$M[1] = initial$M
  mcmc0$Sigma[,,1] = initial$Sigma
  for (i in 1:N0) mcmc0$muh[,,i,1] = initial$mh
  mcmc0$wh[,1] = initial$wh
  mcmc0$r[,1] = initial$r
  lfit1 = survreg(Surv(full_data[id0,1],full_data[id0,3])~cov0,dist="gaussian")
  lfit2 = survreg(Surv(full_data[id0,2],full_data[id0,4])~cov0,dist="gaussian")
  mcmc0$betah1[,,1] = matrix(rep(lfit1$coefficients[-2], H),H,Ncov, byrow=T)
  mcmc0$betah2[,,1] = matrix(rep(lfit2$coefficients[-2], H),H,Ncov, byrow=T)
  mcmc0$imputey[,,1] = full_data[id0,1:2]
  prior_betah1 = mcmc0$betah1[1,,1]
  prior_betah2 = mcmc0$betah2[1,,1]
  
  print("MCMC Iteration under control:")
  for (iter in 2:Niter)
  {
    print(iter)
    tmp = update_r(mcmc0$Sigma[,,iter-1], mcmc0$muh[,,,iter-1],mcmc0$wh[,iter-1], full_data[id0,], mcmc0$r[,iter-1], N0)
    mcmc0$r[,iter] = tmp$r
    mcmc0$imputey[,,iter] = tmp$impute_y
    tmp2 = update_wh_and_M(mcmc0$r[,iter], mcmc0$M[iter-1],lambda1, lambda2)
    mcmc0$M[iter] = tmp2$M
    mcmc0$wh[,iter] = tmp2$wh
    mcmc0$muh[,,,iter] = update_muh(mcmc0$wh[,iter-1],mcmc0$r[,iter], mcmc0$Sigma[,,iter-1], mcmc0$muh[,,,iter-1], mcmc0$imputey[,,iter], mcmc0$betah1[,,iter-1],mcmc0$betah2[,,iter-1], N0, cov0,covariance0, inv_covariance0,partial_cov0)
    mcmc0$betah1[,,iter] = update_betah1(mcmc0$muh[,,,iter], mcmc0$r[,iter],prior_betah1,var0_betah,partial_cov0,Ncov)
    mcmc0$betah2[,,iter] = update_betah2(mcmc0$muh[,,,iter], mcmc0$r[,iter],prior_betah2,var0_betah,partial_cov0,Ncov)
    mcmc0$Sigma[,,iter] = update_Sigma(mcmc0$r[,iter], mcmc0$imputey[,,iter], mcmc0$muh[,,,iter], N0, Phi0,lambda0)
  }
  
  id <- seq(burn.in+1, Niter, lag)
  mcmc1_save <- NULL
  mcmc1_save$M = mcmc1$M[id] 
  mcmc1_save$Sigma = mcmc1$Sigma[,,id] 
  mcmc1_save$muh = mcmc1$muh[,,,id]
  mcmc1_save$betah1 = mcmc1$betah1[,,id]
  mcmc1_save$betah2 = mcmc1$betah2[,,id]
  mcmc1_save$r = mcmc1$r[,id]
  mcmc1_save$wh = mcmc1$wh[,id]
  mcmc1_save$imputey = mcmc1$imputey[,,id]
  
  mcmc0_save <- NULL
  mcmc0_save$M = mcmc0$M[id] 
  mcmc0_save$Sigma = mcmc0$Sigma[,,id] 
  mcmc0_save$muh = mcmc0$muh[,,,id]
  mcmc0_save$betah1 = mcmc0$betah1[,,id]
  mcmc0_save$betah2 = mcmc0$betah2[,,id]
  mcmc0_save$r = mcmc0$r[,id]
  mcmc0_save$wh = mcmc0$wh[,id]
  mcmc0_save$imputey = mcmc0$imputey[,,id]
  
  mcmc0 = mcmc0_save
  mcmc1 = mcmc1_save
  #save(mcmc0, mcmc1, file="saved_mcmc.RData")
  return(list(mcmc0=mcmc0, mcmc1=mcmc1))
}



