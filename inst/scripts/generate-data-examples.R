##
N    <- 6
n    <- 5
beta <- c(-2, 1, 0.5, 0.25)

##
genDataLongV1(N, n, beta, bij=genV0(N, n, sigma=1))
genDataLongV1(N, n, beta, bij=genV1(N, n, sigma=1, lambda=1))
genDataLongV1(N, n, beta, bij=genV2(N, n, sigma0=1, sigma1=1))
genDataLongV1(N, n, beta, bij=genV3(N, n, Sigma=diag(c(1, 1))))
genDataLongV1(N, n, beta, bij=genV4(N, n, sigma=1, rho=0.7))
genDataLongV1(N, n, beta, bij=genV5(N, n, sigma=1))
##
genDataWide(N, n, beta, bij=genV0(N, n, sigma=1))

##
N    <- 6
n    <- 1:6
beta <- c(-2, 1, 0.5, 0.25)
##
genDataLongV2(N, n, beta, bij=0)
