logit <- function(p) {
  log(p/(1-p))
}

expit <- function(x) {
  exp(x)/(1 + exp(x))
}
