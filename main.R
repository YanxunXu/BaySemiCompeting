library(MASS)
library(survival)
library(mc2d)
library(tmvtnorm) #for multivariate normal
library(mvtnorm)
library(MCMCpack) #for inverse whishart
library(mvnfast) #for fast mvrnorm
source("update_pos.R")  #store all MCMC updates functions
source("mcmc.R")  #store the function to run MCMC
source("Marginal_Surv.R")  #store the function to compute the marginal survivals
source("Estimand.R")  #store the function to compute the causal estimand h(u)


##full_data: The first column is observed progression time, the second column is the observed
##survival time, the third column is delta (the censoring indicator for progress), 
## and the last column is the xi (the censoring indicator for survival)
##delta=xi=1, we observe both progression time and death time
##delta=xi=0, we only have the censoring time C
##delta=1, xi=0, we observe progression time, but censored before observing death
##delta=0, xi=1, we observe death only, but neither progression nor death
load("data.Rdata")
########################################################################
Niter = 5000 #the number of MCMC iterations
burn.in = 2000 #the number of burn.in
lag=10 #save samples every lag iterations
nsave = (Niter-burn.in)/lag #the number of saved samples
#Not run. The MCMC result is saved in the file "saved_mcmc.RData"
#mcmc_result <- main_mcmc(Niter, burn.in, lag, full_data, Z)

#Reproduce Figure 1 in Supplementary Material
##############Marginal 
load("saved_mcmc.RData")
#NOT RUN. I save the results below in "Figure1.RData". 
# figure1_result<-Marginal_survival(full_data, Z, mcmc0, mcmc1)
# fmean0_ave=figure1_result$fmean0_ave
# fmean1_ave=figure1_result$fmean1_ave
# fquantile0=figure1_result$fquantile0
# fquantile1=figure1_result$fquantile1

load("Figure1.RData")
library(survival)
pdf("brain_survival.pdf")
tim = seq(0,10,0.3)
mini.surv <- survfit(Surv(full_data[,2], full_data[,4])~ Z[,1], conf.type="none")
plot(mini.surv, col=c(1, 2), xlab="log(Time)", ylab="Survival Probability")
lines(tim, fmean0_ave, col=1,lty=2)
lines(tim, fquantile0[1,], col=1, lty=3)
lines(tim, fquantile0[2,], col=1, lty=3)
lines(tim, fmean1_ave, col=2,lty=2)
lines(tim, fquantile1[1,], col=2, lty=3)
lines(tim, fquantile1[2,], col=2, lty=3)
legend("topright", c("Treatment", "Control"), col = c(2,1),lty=1)
dev.off()
##########################################################
#Reproduce Figure 2 in Supplementary Material
##############h(u)
#NOT RUN. I save hu_result in the file "hu_rho1.RDat" when rho=0.2; and hu_result1 in the file "hu_rho2.RDat" when rho=0.8.  
#rho = 0.2; hu_result <- Estimate_hu(rho, full_data, Z, mcmc0, mcmc1)
#rho = 0.8;hu_result1 <- Estimate_hu(rho, full_data, Z, mcmc0, mcmc1)
load("hu_rho1.RData")
load("hu_rho2.RData")
#####Figure 7 h(u)
u_range = seq(0,6,0.5)
ratio = matrix(0, length(u_range), 2)
ratio_interval=array(0, c(2, length(u_range), 2))
tmp = ifelse(hu_result[,,1]<=0.05 & hu_result[,,2]<=0.05, 1, hu_result[,,2]/hu_result[,,1])
tmp1 = ifelse(hu_result1[,,1]<=0.05 & hu_result1[,,2]<=0.05, 1, hu_result1[,,2]/hu_result1[,,1])

  for (u_index in 1:length(u_range))
  {
    ratio[u_index,1] = mean(tmp[u_index,])
    ratio_interval[, u_index, 1] = quantile(tmp[u_index,], c(0.025, 0.975))
    ratio[u_index,2] = mean(tmp1[u_index,])
    ratio_interval[, u_index, 2] = quantile(tmp1[u_index,], c(0.025, 0.975))
  }

pdf("hu.pdf")
par(mar=c(5.1,5.1,4.1,2.1))
plot(u_range, predict(loess(ratio[,1]~u_range)),"l",xlab="log(Time) (u)", ylab=expression(paste("Estimate of", " ", tau, "(u)")),ylim=c(0, 2),cex.axis=2, cex.lab=2)
lines(u_range, predict(loess(ratio_interval[1,,1]~u_range)),col=1,lty=2)
lines(u_range, predict(loess(ratio_interval[2,,1]~u_range)),col=1,lty=2)
lines(u_range, predict(loess(ratio[,2]~u_range)),col=2)
lines(u_range, predict(loess(ratio_interval[1,,2]~u_range)),col=2,lty=2)
lines(u_range, predict(loess(ratio_interval[2,,2]~u_range)),col=2,lty=2)
legend("topright", c(expression(paste(rho, "=0.2")),expression(paste(rho, "=0.8"))), col = c(1, 2),lty=c(1,1,1))
dev.off()

