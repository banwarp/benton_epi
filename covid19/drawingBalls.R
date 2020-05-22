# drawingBalls.R

# A simple function to model mixing with a small number of infectious persons and daily random interactions.

drawingBalls <- function(n,k,alpha){
  p <- 0
  for(i in 1:alpha) {
    p <- p + choose((n-i),(k-1))
  }
  p <- p/choose(n,k)
  return(p)
}

drawingBallsDistribution <- function(n,k,alpha,N){
  p <- 0
  for(i in 1:alpha) {
    p <- p + choose((n-i),(k-1))
  }
  p <- p/choose(n,k)
  sd <- sqrt(N*p*(1-p))
  return(list(N*p,sd))
}
