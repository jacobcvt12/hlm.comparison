##
#### Potential Scale Reduction
##

##
CalcPSR <- function(ChainMat)
{
  ## number of scans
	n <- dim(ChainMat)[1]
  
	## Between and Within-chain variation
	B <- n * var(apply(ChainMat, 2, mean))
  W <- mean(apply(ChainMat, 2, var))

  ## MPV = marginal posterior variance
	MPV <- ((n-1)*W + B) / n

  ## psr = potential scale reduction
  PSR <- sqrt(MPV / W)
  return(PSR)
}

