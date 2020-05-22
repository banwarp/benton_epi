# drawingBalls.R

drawingBalls <- function(n,k,alpha){
  p <- 0
  for(i in 1:alpha) {
    p <- p + choose((n-i),(k-1))
  }
  p <- p/choose(n,k)
  return(p)
}