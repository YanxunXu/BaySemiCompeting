#response = data[,1:2]
init <- function(response, Num_Patient, H)
{
  hc = hclust(dist(response)^2, "cen")
	r = cutree(hc, k=H)
  wh1 <- table(r)/Num_Patient
  idx <- order(wh1,decreasing=T)
  wh <- wh1[idx]
  tmp0 = split(response,r)
  mh = matrix(NA, H, 2)
  for (h in 1:H)
  {
    mh[h,] = colSums(matrix(tmp0[[h]],ncol=2))/sum(r==h)
  }
  M = 1
  Sigma = matrix(c(0.5, 0.1, 0.1, 0.5),2,2)
  return(list(mh=mh,wh=wh, M=M,Sigma=Sigma, r=r))
}

#wh=mcmc0$wh[,iter-1]; r=mcmc0$r[,iter]; muh=mcmc0$muh[,,,iter-1];data=full_data[id0,];Sigma=mcmc0$Sigma[,,,iter-1]; 
#cov=cov0; covariance=covariance0; inv_covariance=inv_covariance0; r=mcmc0$r[,iter-1];Num_Patient=N0;
update_r <- function(Sigma, muh, wh, data, r, Num_Patient)
{
  impute_y = data[,1:2]
	for (i in 1:Num_Patient)
	{
    if (data[i,3]==1 & data[i,4]==1)
		{
      tmp0 = rep(0, H)
			ph = diag(exp(-0.5*t(data[i,1:2]-t(muh[,,i]))%*%solve(Sigma)%*%(data[i,1:2]-t(muh[,,i]))))*wh
        #dmvnorm(data[i,1:2], muh[2,], Sigma)*wh
		  r[i] = sample(1:H, 1, prob=ph)
		}
		if (data[i,3]==0 & data[i,4]==0)
		{
      tmp0 = rep(0, H)
      for (h in 1:H) tmp0[h] = pmvnorm(lower=data[i,1:2], upper=c(Inf,Inf),mean=muh[h,,i], sigma=Sigma)
			ph = tmp0*wh
      ph[ph<=0] = 0
		  r[i] = sample(1:H, 1, prob=ph)
			impute_y[i,] = rtmvnorm(1, mean=muh[r[i],,i], sigma=Sigma, lower=data[i,1:2],algorithm="gibbs")
			if (impute_y[i,1]==Inf | impute_y[i,2]==Inf) impute_y[i,]=c(data[i,1],data[i,2])
		}
		if (data[i,3]==1 & data[i,4]==0)
		{
		  tmp0 = rep(0, H)
		  for (h in 1:H) tmp0[h] = dtmvnorm.marginal(data[i,1],n=1, mean=muh[h,,i],sigma=Sigma, lower=c(-Inf,data[i,2]),upper=c(Inf,Inf))
		  tmp0[tmp0=="NaN"] = 0
		  tmp0[tmp0=="Inf"] = 0
      ph = tmp0*wh
		  r[i] = sample(1:H, 1, prob=ph)
		  impute_y[i,] = rtmvnorm(1, mean=muh[r[i],,i], sigma=Sigma, lower=c(data[i,1],data[i,2]),upper=c(data[i,1]+1e-6,Inf),algorithm="gibbs")
		  if (impute_y[i,1]==Inf | impute_y[i,2]==Inf) impute_y[i,]=c(data[i,1],data[i,2])
		}
		if (data[i,3]==0 & data[i,4]==1)
		{
		  tmp0 = rep(0, H)
		  for (h in 1:H) tmp0[h] = dtmvnorm.marginal(data[i,2],n=2, mean=muh[h,,i],sigma=Sigma, lower=c(data[i,1],-Inf),upper=c(Inf,Inf))
		  tmp0[tmp0=="NaN"] = 0
		  tmp0[tmp0=="Inf"] = 0
      ph = tmp0*wh
		  r[i] = sample(1:H, 1, prob=ph)
		  impute_y[i,] = rtmvnorm(1, mean=muh[r[i],,i], sigma=Sigma, lower=c(data[i,1],data[i,2]),upper=c(Inf,data[i,2]+1e-6),algorithm="gibbs")
		  if (impute_y[i,1]==Inf | impute_y[i,2]==Inf) impute_y[i,]=c(data[i,1],data[i,2])
		}
	}
	return(list(r=r,impute_y=impute_y))
}

#r<-mcmc$r[,iter]; M=mcmc$M[iter-1]
update_wh_and_M <- function(r, M, lambda1, lambda2)
{
   ## returns: wh
    vh <- rep(0,H)  # initialize
    wh <- rep(0,H)
    V <-  1         # record prod_{g<h} (1-vh_h)
    for(h in 1:(H-1)){
      Ah <- which(r==h)
      Bh <- which(r>h)
      vh[h] <-  rbeta(1, 1+length(Ah), M+length(Bh))
      wh[h] <- vh[h]*V
      V <- V*(1-vh[h])
    }
    vh[H] <- 1.0
    wh[H] <- V
    #M <- rgamma(1, lambda1+H-1, lambda2-sum(log(1-vh[1:(H-1)])) )
    M <- rgamma(1, lambda1+H-1, lambda2-max(-100, sum(log(1-vh[1:(H-1)]))) )
    #log(.Machine$double.xmin) = -708.3964
    return(list(wh=wh, M=M))
}

#wh=mcmc0$wh[,iter-1]; r=mcmc0$r[,iter]; muh=mcmc0$muh[,,,iter-1];response=mcmc0$imputey[,,iter];Sigma=mcmc0$Sigma[,,iter-1]; betah1=mcmc0$betah1[,,iter-1];betah2=mcmc0$betah2[,,iter-1]
#cov=cov0; covariance=covariance0; inv_covariance=inv_covariance0; r=mcmc0$r[,iter-1];Num_Patient=N0;partial_cov=partial_cov0
update_muh <- function(wh,r,Sigma, muh, response, betah1,betah2, Num_Patient, cov, covariance, inv_covariance,partial_cov)
{
	   for(h in 1:H){
	     if(any(r==h)){      # some data assigned to h-th pointmass
	       Sh <- which(r==h) 
	       nh <- length(Sh)
         W = matrix(0, nh, Num_Patient)
         for (i in 1:nh)
         {
           W[i, Sh[i]] = 1
         }
         response1 = response[Sh,1] - Sigma[1,2]/Sigma[2,2]*(response[Sh,2]-muh[h, 2, Sh])
         sigma21 = Sigma[1,1] - Sigma[1,2]*Sigma[2,1]/Sigma[2,2]
	       var1 = chol2inv(chol(inv_covariance+1/sigma21*t(W)%*%W))
	       mu1 = var1%*%(t(W)%*%response1/sigma21 + t(partial_cov)%*%betah1[h,])
	       muh[h,1,] <- rmvn(1, t(mu1), var1)
	       response2 = response[Sh,2] - Sigma[2,1]/Sigma[1,1]*(response[Sh,1]-muh[h, 1, Sh])
	       sigma22 = Sigma[2,2] - Sigma[1,2]*Sigma[2,1]/Sigma[1,1]
	       var2 = chol2inv(chol(inv_covariance+1/sigma22*t(W)%*%W))
	       mu2 = var2%*%(t(W)%*%response2/sigma22 + t(partial_cov)%*%betah2[h,])
	       muh[h,2,] <- rmvn(1, t(mu2), var2)
	     } else {            # no data assinged to h-th pointmass# sample from base measure
	       mu1 <- cov%*%betah1[h,]
	       mu2 <- cov%*%betah2[h,]
	       muh[h,1,] <- rmvn(1, t(mu1), covariance)
	       muh[h,2,] <- rmvn(1, t(mu2), covariance)
	     }
	   }
    return(muh)
}


#muh <- mcmc1$muh[,,,iter];r=mcmc1$r[,iter] 
update_betah1 <- function(muh, r,prior_betah1, var,partial_cov, Ncov)
{
	betah1 <- matrix(0,H, Ncov)     # initialize
     for(h in 1:H){
      if(any(r==h)){      # some data assigned to h-th pointmass
        mu = var%*%(partial_cov%*%muh[h,1,] + diag(Ncov)%*%prior_betah1 )
        betah1[h,] <- rmvn(1, t(mu), var)
      } else {            # no data assinged to h-th pointmass     # sample from base measure
        betah1[h,] <- rmvn(1, prior_betah1, diag(Ncov))}
    }
    return(betah1)
}





#muh <- mcmc$muh[,,,iter];r=mcmc$r[,iter] 
update_betah2 <- function(muh, r,prior_betah2, var, partial_cov, Ncov)
{
  betah2 <- matrix(0,H, Ncov)     # initialize
  for(h in 1:H){
    if(any(r==h)){      # some data assigned to h-th pointmass
      mu = var%*%(partial_cov%*%muh[h,2,] + diag(Ncov)%*%prior_betah2 )
      betah2[h,] <- rmvn(1, t(mu), var)
    } else {            # no data assinged to h-th pointmass     # sample from base measure
      betah2[h,] <- rmvn(1, prior_betah2, diag(Ncov))}
  }
  return(betah2)
}


#mcmc=mcmc0; muh=mcmc$muh[,,,iter];r=mcmc$r[,iter]; response=mcmc$imputey[,,iter];Num_Patient=N0; Phi=Phi0
update_Sigma <- function(r, response, muh, Num_Patient, Phi,lambda0)
{
	tmp = matrix(0, Num_Patient, 2)
	for (h in 1:H)
	{
		if(any(r==h))
		{
			Sh = which(r==h)
			tmp[Sh,] = t(muh[h,,Sh])
		}
	}
	Sigma <- riwish(lambda0+Num_Patient, Phi + t(response-tmp)%*%(response-tmp))
	return(Sigma)
}

#wh<-mcmc$wh[,iter]; muh=mcmc$muh[,,iter]; sigma2=mcmc$sigma2[iter];r=mcmc$r[,iter];covD=cov;
fmean <- function(wh, muh, theta, cov, covD, inv_covariance, Npat, sigma2)
{
  fx = rep(0, Npat)
  for (i in 1:Npat)
  {
    coeff = t(covD)-cov[i,]
    c1 = exp(-colSums(coeff^2*theta^2))
    fx[i] = exp(sum(wh*c1%*%inv_covariance%*%t(muh)+1/2*wh^2*sigma2))
  }
  return(fx)
}




#xgrid= tim;wh=mcmc0$wh[,iter]; muh=mcmc0$muh[,2,,iter]; betah=mcmc0$betah2[,,iter];sigma2=mcmc0$Sigma[2,2,iter];dpat=covD
#xgrid= tim;wh=mcmc1$wh[,iter]; muh=mcmc1$muh[,1,,iter]; betah=mcmc1$betah1[,,iter];sigma2=mcmc1$Sigma[1,,iter]
#wh=mcmc1$wh[,iter]; muh=mcmc1$muh[,2,,iter]; betah=mcmc1$betah2[,,iter];sigma2=mcmc1$Sigma[2,2,iter]
fmar_survival <- function(xgrid, wh, muh, betah, sigma2, dpat,cov,inv_covariance)
{
  fx <- rep(0, length(xgrid))
  for (h in 1:H)
  {
    coeff = t(cov)-dpat
    c1 = exp(-colSums(coeff^2))
    var1 = 1.01 - c1%*%inv_covariance%*%c1
    fx <- fx + wh[h]*pnorm(xgrid, m=betah[h,]%*%dpat+c1%*%inv_covariance%*%(muh[h,]-cov%*%betah[h,]), sd=sqrt(sigma2 + var1)) 
  }
  return(fx)
}







fbar_iter <- function(mcmc0, mcmc1, u, rho, Npat,partial_cov0,partial_cov1, cov0, cov1, cov,var1, var0)
{
  Mrep = 100
  hu = array(0, c(nsave,Npat,2))
  hu_ave = matrix(0, nsave, 2)
  for (iter in 1:nsave)
  {
    wh0=mcmc0$wh[,iter]; muh0=mcmc0$muh[,,,iter];
    betah10=mcmc0$betah1[,,iter];betah20=mcmc0$betah2[,,iter]; 
    Sigma0=mcmc0$Sigma[,,iter];
    wh1=mcmc1$wh[,iter]; muh1=mcmc1$muh[,,,iter];
    betah11=mcmc1$betah1[,,iter];betah21=mcmc1$betah2[,,iter]; 
    Sigma1=mcmc1$Sigma[,,iter]; 
     
    mu10_all = betah10%*%t(cov) + t(partial_cov0%*%(t(muh0[,1,])-cov0%*%t(betah10)) )
    mu20_all = betah20%*%t(cov) + t(partial_cov0%*%(t(muh0[,2,])-cov0%*%t(betah20)) )
    
    mu11_all = betah11%*%t(cov) + t(partial_cov1%*%(t(muh1[,1,])-cov1%*%t(betah11)) )
    mu21_all = betah21%*%t(cov) + t(partial_cov1%*%(t(muh1[,2,])-cov1%*%t(betah21)) )
    
    for (k in 1:Npat)
    {
      tmp1 = sample(1:H, Mrep, replace=T, prob=wh1)
      mean_tmp1 = t(rbind(mu11_all[,k], mu21_all[,k]))
      mean1 = mean_tmp1[tmp1,]
      tem_sample1 = matrix(0, Mrep, 2)
      tem_sample1 = mgcv::rmvn(Mrep, mean1, Sigma1 + var1[k]*diag(2))
      temp2 = 0
      for (h in 1:H)
      {
        temp2 = temp2 + wh1[h]*pnorm(tem_sample1[,2], mu21_all[h,k], sqrt(Sigma1[2,2] + var1[k]))
      }
      tmp2 = rho*qnorm(temp2) + rnorm(Mrep, 0, sqrt(1-rho^2))
      tmp_value = sum(wh0*pnorm(u, mu20_all[,k], sqrt(Sigma0[2,2] + var0[k])) )
      if (tmp_value>0.999) tmp_value=0.999
      tmp3 = qnorm( tmp_value )
      num1 = sum(tem_sample1[,1] < u & tem_sample1[,2] >=u & tmp2 >= tmp3)
      
      tmp0 = sample(1:H, Mrep, replace=T, prob=wh0)
      mean_tmp0 = t(rbind(mu10_all[,k], mu20_all[,k]))
      mean0 = mean_tmp0[tmp0,]
      tem_sample0 = matrix(0, Mrep, 2)
      tem_sample0 = mgcv::rmvn(Mrep, mean0, Sigma0 + var0[k]*diag(2))
      temp4 = 0
      for (h in 1:H)
      {
        temp4 = temp4 + wh0[h]*pnorm(tem_sample0[,2], mu20_all[h,k], sqrt(Sigma0[2,2] + var0[k]))
      }
      tmp4 = rho*qnorm(temp4)
      tmp4 = tmp4 + rnorm(Mrep, 0, sqrt(1-rho^2))
      tmp_value = sum(wh1*pnorm(u, mu21_all[,k], sqrt(Sigma1[2,2] + var1[k])) )
      if (tmp_value>0.999) tmp_value=0.999
      tmp5 = qnorm( tmp_value )
      num0 = sum(tem_sample0[,1] < u & tem_sample0[,2] >=u & tmp4 >= tmp5)
      hu[iter,k,1] = num0
      hu[iter,k,2] = num1
    }
  }
    hu_ave[,1] = apply(hu[,,1], 1, mean)
    hu_ave[,2] = apply(hu[,,2], 1, mean)
  return(hu_ave)
}

