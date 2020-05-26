# drawingBalls.R

# Note: This script does not use the built-in choose(n,k) but instead computes n choose k from the sum/difference of logs.
# This is a bit slower for large n and k, but it can admit much larger n and k then choose(n,k).

drawingBalls <- function(n,k,alpha){
  p <- 0
  j <- k-1
  for(i in 1:alpha) {
    m <- n-i
    top <- c((m-j+1):m)
    sumlogtop <- sum(log(top))
    bottom <- c(1:j)
    sumlogbottom <- sum(log(bottom))
    padd <- exp(sumlogtop-sumlogbottom)
    p <- p + padd
  }
  top <- c((n-k+1):n)
  sumlogtop <- sum(log(top))
  bottom <- c(1:k)
  sumlogbottom <- sum(log(bottom))
  p <- log(p)-(sumlogtop-sumlogbottom)
  p <- exp(p)
  return(p)
}

drawingBallsDistribution <- function(n,k,alpha,N){
  p <- drawingBalls(n,k,alpha)
  sd <- sqrt(N*p*(1-p))
  return(list(mean = N*p,stdev = sd))
}

drawingBallsSamples <- function(n,k,alpha,N,nTrials) {
  p <- drawingBalls(n,k,alpha)
  return(rbinom(nTrials,N,p))
}
