#' Generates random effects data
#'
#' Function can generate data with random intercept and fixed slope,
#' or random intercept and slope with arbitrary fixed covariates
#' @param clusters The number of clusters (> 1)
#' @param cluster.size Either a vector indicating the size of each cluster or a single number which indicates the size of each cluster (same across clusters)
#' @param params.fe A list containing the values of the fixed effects
#' @param params.re A list contiaining the values of the random effects
#' @param G "True" mixing distribution as a character (e.g. "normal", "gamma")
#' @param ... Named parameters for the mixing distribution G
#'
#' @return data.frame
#'
#' @export
create.data <- function(clusters,
                        cluster.size,
                        params.fe,
                        params.re,
                        G,
                        ...) {
  return(NULL)
}

##
## 0	Normally distributed random intercepts with common standard deviation
## 1	Non-Normal random interceps; gamma(lambda, 1) distribution
## 2	Random intercepts with standard deviation depending on value of a (binary) cluster-level covariate
## 3	Random intercepts and slopes
## 4	Autocorrelated random effects
## 5  Two-point distribution
##

##
genV0 <- function(N, n, sigma)
{
	## Normally distributed random intercepts with common standard deviation
	##
	bi  <- rnorm(N, mean=0, sd=sigma)
	bij <- rep(bi, rep(n, N))

	return(bij)
}

##
genV1 <- function(N, n, sigma, lambda)
{
	## Non-Normal random interceps; gamma(lambda, 1) distribution
	##
	ai  <- rgamma(N, lambda, 1)
	bi  <- sigma * (ai - lambda) / sqrt(lambda)
	bij <- rep(bi, rep(n, N))

	return(bij)
}

##
genV2 <- function(N, n, sigma0, sigma1)
{
	## Random intercepts with standard deviation depending on value of a (binary) cluster-level covariate
	##
	bi0 <- rnorm(N/2, mean=0, sd=sigma0)
	bi1 <- rnorm(N/2, mean=0, sd=sigma1)
	bij <- c(rep(bi0, rep(n, N/2)), rep(bi1, rep(n, N/2)))

	return(bij)
}

##
genV3 <- function(N, n, Sigma)
{
	## Observation-specific covariate
	#Xij2 <- rep((c(1:n)-1)/(n-1), N)
	#Xij2 <- rnorm(n*N, 0, 2)
	Xij2 <- rep(c(1:n)-median(1:n), N)

	## Random intercepts and slopes
	##
	bi  <- rmvnorm(N, mean=c(0,0), sigma=Sigma)
	bij <- rep(bi[,1], rep(n, N)) + (Xij2 * rep(bi[,2], rep(n, N)))

	return(bij)
}

##
genV4 <- function(N, n, sigma, rho)
{
	## Autocorrelated random effects
	##
	covStruct <- matrix(sigma^2, nrow=n, ncol=n)
	for(j in 1:n)
	{
		for(k in 1:n) covStruct[j,k] <- covStruct[j,k] * (rho^(abs(j-k)))
	}

	##
	bij <- c(t(rmvnorm(N, mean=rep(0,n), sigma=covStruct)))

	return(bij)
}

##
genV5 <- function(N, n, sigma)
{
	## Two-point, symmetric distribution with variance sigma^2
	##
	bVal <- sigma/sqrt(2)
	bi   <- sample(c(-bVal, bVal), N, replace=TRUE)
	bij  <- rep(bi, rep(n, N))

	return(bij)
}

##
genDataLongV1 <- function(N, n, beta, bij)
{
	## Version V1 holds 'n' fixed across clusters

	## Cluster-specific covariate
	Xij1 <- rep(c(0,1), rep((N*n)/2, 2))

	## Observation-specific covariate
	# Xij2 <- rep((c(1:n)-1)/(n-1), N)
	# Xij2 <- rnorm(n*N, 0, 2)
	Xij2 <- rep(c(1:n)-median(1:n), N)

	## Conditionally specific linear predictor: \mu_{ij}^c
	muC <- beta[1] + (beta[2] * Xij1) + (beta[3] * Xij2) + bij
	Yij <- rbinom(N*n, 1, expit(muC))

	##
	value <- as.data.frame(cbind(rep(1:N, rep(n, N)), Yij, Xij1, Xij2, bij))
  names(value) <- c("id", "Y", "X1", "X2", "b")
	return(value)
}

##
genDataLongV2 <- function(N, n, beta, bij)
{
	## Version V2 permits 'n' to vary across clusters

	## Cluster-specific covariate
	Xij1 <- rep(c(0,1), c(sum(n[c(1:(N/2))]), sum(n[-c(1:(N/2))])))

	## Observation-specific covariate
	Xij2 <- rnorm(sum(n), 0, 2)

	## Conditionally specific linear predictor: \mu_{ij}^c
	muC <- beta[1] + (beta[2] * Xij1) + (beta[3] * Xij2) + (beta[4] * Xij1 * Xij2) + bij
	Yij <- rbinom(sum(n), 1, expit(muC))

	##
	value <- as.data.frame(cbind(rep(1:N, n), Yij, Xij1, Xij2, Xij1*Xij2, bij))
  names(value) <- c("id", "Y", "X1", "X2", "X12", "b")
	return(value)
}

##
genDataWide <- function(N, n, beta, bij)
{
	## Cluster-specific covariate
	X1   <- rep(c(0,1), rep(N/2, 2))
	Xij1 <- rep(c(0,1), rep((N*n)/2, 2))

	## Observation-specific covariate
	X2   <- c(1:n - 1) / (n-1)
	Xij2 <- rep(c(1:n - 1)/(n-1), N)

	## Conditionally specific linear predictor: \mu_{ij}^c
	muC <- beta[1] + (beta[2] * Xij1) + (beta[3] * Xij2) + (beta[4] * Xij1 * Xij2) + bij
	Yij <- matrix(rbinom(N*n, 1, expit(muC)), nrow=N, byrow=TRUE)

	##
	value <- list(id=c(1:N), Yij=Yij, X1=X1, X2=X2, b=bij[seq(from=1, to=(N*n), by=n)])
	return(value)
}
