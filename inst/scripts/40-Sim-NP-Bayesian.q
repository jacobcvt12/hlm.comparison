
## Includes a total of 9 analyses (7 of which are Bayesian)
##   - naive frequentist analysis via glm()
##   - frequentist analysis via glmmML()
##   - Bayesian analysis via LogisticNormal()
##   - Bayesian analysis via DPglmm()
##       * 6 of these
##
## Random effect values taken as the quantiles from the appropriate distribution
##  - remain the same across simulated datasets
##  - use those that were generate in Simulation #2
##
## Only consider N=100 and n=10


##
rm(list=ls())

##
index <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#index <- as.numeric(Sys.getenv('LSB_JOBINDEX'))
##
library(glmmML, lib="~/Rpackages")
library(DPpackage, lib="~/Rpackages")
##
R <- 10

##
seedBase <- 1234
set.seed(seedBase + index)

##

#source("/Users/josephantonelli/Documents/transfer_files/GLMMs/Utilities/SimulateData.q")
#source("/Users/josephantonelli/Documents/transfer_files/GLMMs/Utilities/PSR.q")
#source("/Users/josephantonelli/Documents/transfer_files/GLMMs/Bayes/LogisticNormal/LogisticNormal.q")

source("~/GLMMs/Utilities/SimulateData.q")
source("~/GLMMs/Utilities/PSR.q")
source("~/GLMMs/Bayes/LogisticNormal/LogisticNormal.q")
expit <- function(x) exp(x)/(1 + exp(x))

gamm = function(N2,lambda, sigma) {
  gam = rgamma(N2,lambda,1)
  random = sigma * (gam - lambda) / sqrt(lambda)  # random intercepts
  return(random)
}

##
N     <- 200
n     <- 10
betaV <- c(-2, 1, 0.5)
##

##
nChains <- 3
nScans  <- 10000
thin    <- 4
burnin  <- (nScans/thin) * 0.2
nsave   <- (nScans/thin) - burnin
mcmc    <- list(nburn=burnin, nsave=nsave, nskip=thin, ndisplay=500)
##
estFE <- array(NA, dim=c(R,7,5,14))
estRE <- array(NA, dim=c(R,N,5,15))   ## final element of the 4th dimension stores the actual values of the REs
##
psrFE <- array(NA, dim=c(R,7,5,14))  
psrRE <- array(NA, dim=c(R,N,5,14))

seFE <- array(NA, dim=c(R,3,5,14))
EBconverge <- array(1, dim=c(R,5))
 
##
for(r in 1:R)
{
	trueRE <- matrix(NA, nrow=N, ncol=5)
	trueRE[,1] <- rnorm(N, mean=0, sd=2)
	trueRE[,2] <- gamm(N, 0.5, 2)
	trueRE[,3] <- c(rnorm(N/2, mean=0, sd=sqrt(7)), rnorm(N/2, mean=0, sd=1))
	trueRE[,4] <- sample(c(-2,2), N, replace=TRUE)
	trueRE[,5] <- rt(N, df=8/3)

	##
	for(type in 1:5)
	{	
		##
		estRE[r, ,type,dim(estRE)[4]] <- trueRE[,type]
		if(type == 1) bij <- rep(trueRE[,1], rep(n, N))
		if(type == 2) bij <- rep(trueRE[,2], rep(n, N))
		if(type == 3) bij <- rep(trueRE[,3], rep(n, N))
		if(type == 4) bij <- rep(trueRE[,4], rep(n, N))
		if(type == 5) bij <- rep(trueRE[,5], rep(n, N))

		##
		simData  <- genDataLongV1(N, n, betaV, bij)
		
		## Naive frequentist analysis
		##
		Aindex <- 1
		##
		fitNaive                 <- glm(Y ~ X1 + X2, family=binomial, data=simData)
		estFE[r,1:3,type,Aindex] <- fitNaive$coefficients
		seFE[r,,type,Aindex] <- sqrt(as.numeric(diag(vcov(fitNaive))))

		## Frequentist Logistic-Normal analysis
		##
		Aindex <- 2
		##
		freqFit                  <- glmmML(Y ~ X1 + X2, cluster=id, data=simData, family=binomial)
		estFE[r,1:3,type,Aindex] <- freqFit$coefficients
		estRE[r,,type,Aindex]    <- freqFit$posterior.modes
		seFE[r,,type,Aindex]     <- sqrt(diag(freqFit$variance))[1:3]

		## Bayesian Logistic-Normal analysis
		##
		Aindex <- 3
		##
		LNBayesFE <- array(NA, dim=c(nScans/thin, 4, nChains))
		LNBayesRE <- array(NA, dim=c(nScans/thin, N, nChains))
		##
		for(chain in 1:nChains)
		{
			temp <- LogisticNormal(fitNaive, simData$id, nScans=nScans, thin=thin, startValues=c(fitNaive$coef, 1) * runif(5, 0.9, 1.1), startV=rep(0, N))
			LNBayesFE[,,chain] <- temp[,c(1:4)]
			LNBayesRE[,,chain] <- temp[,-c(1:4)]
		}
		##
		estFE[r,c(1:3,5),type,Aindex] <- apply(LNBayesFE[-c(1:burnin),,], 2, FUN=median)
		estRE[r,,type,Aindex]         <- apply(LNBayesRE[-c(1:burnin),,], 2, FUN=median)
		##
		psrFE[r,c(1:3,5),type,Aindex] <- apply(LNBayesFE[-c(1:burnin),,], 2, FUN=CalcPSR)
		psrRE[r,,type,Aindex]         <- apply(LNBayesRE[-c(1:burnin),,], 2, FUN=CalcPSR)
		
		seFE[r,,type,Aindex]          <- apply(LNBayesFE[-c(1:burnin),,], 2, FUN=sd)[1:3]

		## Bayesian: Dirichlet process
		##
		for(Aindex in 4:14) {
			##
			if(Aindex == 4) prior <- list(beta0=rep(0,2), Sbeta0=diag(1000,2), a0=NULL, b0=0, alpha=100, mub=rep(0,1), Sb=diag(1000,1), nu0=6, tinv=diag(6,1))
			if(Aindex == 5) prior <- list(beta0=rep(0,2), Sbeta0=diag(1000,2), a0=NULL, b0=0, alpha=50, mub=rep(0,1), Sb=diag(1000,1), nu0=6, tinv=diag(6,1))
			if(Aindex == 6) prior <- list(beta0=rep(0,2), Sbeta0=diag(1000,2), a0=NULL, b0=0, alpha=25, mub=rep(0,1), Sb=diag(1000,1), nu0=6, tinv=diag(6,1))
			if(Aindex == 7) prior <- list(beta0=rep(0,2), Sbeta0=diag(1000,2), a0=NULL, b0=0, alpha=10, mub=rep(0,1), Sb=diag(1000,1), nu0=6, tinv=diag(6,1))
			if(Aindex == 8) prior <- list(beta0=rep(0,2), Sbeta0=diag(1000,2), a0=NULL, b0=0, alpha=5, mub=rep(0,1), Sb=diag(1000,1), nu0=6, tinv=diag(6,1))
			if(Aindex == 9) prior <- list(beta0=rep(0,2), Sbeta0=diag(1000,2), a0=NULL, b0=0, alpha=1, mub=rep(0,1), Sb=diag(1000,1), nu0=6, tinv=diag(6,1))
			if(Aindex == 10) prior <- list(beta0=rep(0,2), Sbeta0=diag(1000,2), a0=NULL, b0=0, alpha=0.1, mub=rep(0,1), Sb=diag(1000,1), nu0=6, tinv=diag(6,1))
			## using flat prior centered at N/2
			if(Aindex == 11) prior <- list(beta0=rep(0,2), Sbeta0=diag(1000,2), a0=1.5, b0=0.0125, mub=rep(0,1), Sb=diag(1000,1), nu0=6, tinv=diag(6,1))
			## using kullback liebler with flat prior info on K
			if(Aindex == 12)  prior <- list(beta0=rep(0,2), Sbeta0=diag(1000,2), a0=.491, b0=.004, mub=rep(0,1), Sb=diag(1000,1), nu0=6, tinv=diag(6,1))
			##
			if (Aindex == 13) {
				alpha = 15
				count = 1
				loop=TRUE
				while(loop) {
					nburn <- 0
					nsave <- 1000
					nskip <- 0
					ndisplay <- 500
					mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)
  					prior <- list(beta0=rep(0,2), Sbeta0=diag(1000,2), a0=NULL, b0=0, alpha=alpha[length(alpha)], mub=rep(0,1), 									Sb=diag(1000,1), nu0=6, tinv=diag(6,1))
  
  					fit2 <- try(DPglmm(fixed=Y~X1+X2, random=~1|id, family=binomial(logit), 
  					                   prior=prior, mcmc=mcmc, state=state, status=TRUE, data=simData)) 
  					                    
  					root.alpha = function(x, Ek) {
    					vec.i = 1:N
    					return(Ek - sum(x / (x + vec.i - 1)))
  					}
  
  					if (class(fit2) == "DPglmm") {
  					
  						Ek = mean(fit2$save.state$thetasave[,6])
  
  						temp.alpha = try(uniroot(root.alpha, interval=c(0.00001, 500), Ek=Ek))
  						if (class(temp.alpha) != "list") break
  						new.alpha = temp.alpha$root
  						err = abs(new.alpha - alpha[length(alpha)])
  						alpha = c(alpha,new.alpha)
  						
  						if (length(alpha) > 100) EBconverge[r,type] = 0
  						
  						if (length(alpha) > 2) {
  						    diff1 = sign(alpha[length(alpha)] - alpha[length(alpha) - 1])
    						diff2 = sign(alpha[length(alpha)-1] - alpha[length(alpha) - 2])
    						if (err < 2e-2 | length(alpha) > 50 | diff1 != diff2) break
  						}
 					}
  					if (count > 100)  {
  						EBconverge[r,type] = 0
  						break
  					}
  					count = count + 1
				}
				nChains <- 3
				nScans  <- 10000
				thin    <- 4
				burnin  <- (nScans/thin) * 0.2
				nsave   <- (nScans/thin) - burnin
				mcmc    <- list(nburn=burnin, nsave=nsave, nskip=thin, ndisplay=500)
				prior <- list(beta0=rep(0,2), Sbeta0=diag(1000,2), a0=NULL, b0=0, alpha=alpha[length(alpha)], mub=rep(0,1), 									Sb=diag(1000,1), nu0=6, tinv=diag(6,1))
			}
			
			if (Aindex == 14) {
				aIS = 0.5
				bIS = 0.01
				prior <- list(beta0=rep(0,2), Sbeta0=diag(1000,2), a0=aIS, b0=bIS, mub=rep(0,1), Sb=diag(1000,1), nu0=6, tinv=diag(6,1))
				fit2 <- try(DPglmm(fixed=Y~X1+X2, random=~1|id, family=binomial(logit), 
  					                   prior=prior, mcmc=mcmc, state=state, status=TRUE, data=simData))
  					                   
  				density.gamma = density(fit2$save.state$thetasave[,7], n=8192, from=0.1, to=3*N)

  				weights = 1 / dgamma(density.gamma$x, aIS, bIS)
  				weights_norm = weights / sum(weights)
  
  				density.uniform = density.gamma$y*weights_norm
  				
  			  	m = max(density.uniform)
  				w = which(density.uniform == m)[1]
  				alphaStar = density.gamma$x[w]
  				
  				prior <- list(beta0=rep(0,2), Sbeta0=diag(1000,2), a0=NULL, b0=0, alpha=alphaStar, mub=rep(0,1), Sb=diag(1000,1), nu0=6, tinv=diag(6,1))
			}
			
			##
			DPBayesFE <- array(NA, dim=c(nsave, 7, nChains))
			DPBayesRE <- array(NA, dim=c(nsave, N, nChains))
			##
			chain <- 1
			while(chain < (nChains+1))
			{
				temp <- try(DPglmm(fixed=Y~X1+X2, random=~1|id, family=binomial(logit), prior=prior, mcmc=mcmc, state=state, status=TRUE, data=simData))
				if(class(temp) == "DPglmm")
				{
					DPfix  <- temp$save.state$thetasave
					DPran  <- temp$save.state$randsave
					DPBayesFE[,,chain] <- DPfix
					DPBayesRE[,,chain] <- DPran[,1:N] - matrix(DPfix[,1], nrow=nrow(DPran[,1:N]), ncol=ncol(DPran[,1:N]))
					chain <- chain + 1
				}
			}
			##
			estFE[r,,type,Aindex] <- apply(DPBayesFE, 2, FUN=median)
			estRE[r,,type,Aindex] <- apply(DPBayesRE, 2, FUN=median)
			##
			psrFE[r,,type,Aindex] <- apply(DPBayesFE, 2, FUN=CalcPSR)
			psrRE[r,,type,Aindex] <- apply(DPBayesRE, 2, FUN=CalcPSR)
			
			seFE[r,,type,Aindex]  <- apply(DPBayesFE, 2, FUN=sd)[1:3]
		}

  	##
  	if(type == 1) print(paste(R, r, "Normal done"))
  	if(type == 2) print(paste(R, r, "Standardized Gamma done"))
  	if(type == 3) print(paste(R, r, "Mixture of Normals done"))
  	if(type == 4) print(paste(R, r, "Two-point distribution done"))
   	if(type == 5) print(paste(R, r, "T distribution done"))
	}

	##
	save("estFE", file=paste("~/GLMMs/Misspecification/1_Fixed_N/Output/40_Sim/4-estFE-n", N, "-", n, "-", index, ".dat", sep=""))
	save("estRE", file=paste("~/GLMMs/Misspecification/1_Fixed_N/Output/40_Sim/4-estRE-n", N, "-", n, "-", index, ".dat", sep=""))
	save("psrFE", file=paste("~/GLMMs/Misspecification/1_Fixed_N/Output/40_Sim/4-psrFE-n", N, "-", n, "-", index, ".dat", sep=""))
	save("psrRE", file=paste("~/GLMMs/Misspecification/1_Fixed_N/Output/40_Sim/4-psrRE-n", N, "-", n, "-", index, ".dat", sep=""))
	save("seFE", file=paste("~/GLMMs/Misspecification/1_Fixed_N/Output/40_Sim/4-seFE-n", N, "-", n, "-", index, ".dat", sep=""))
	save("EBconverge", file=paste("~/GLMMs/Misspecification/1_Fixed_N/Output/40_Sim/4-EB-n", N, "-", n, "-", index, ".dat", sep=""))

}
