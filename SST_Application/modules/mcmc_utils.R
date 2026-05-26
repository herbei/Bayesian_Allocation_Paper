rinvgamma1 <- function(shape, rate) {
  1 / rgamma(1, shape = shape, rate = rate)
}

inv_logit <- function(z) {
  if (is.na(z)) return(0.5)
  if (is.infinite(z)) return(if (z > 0) 1 else 0)
  if (z >= 0) {
    1 / (1 + exp(-z))
  } else {
    ez <- exp(z)
    ez / (1 + ez)
  }
}
