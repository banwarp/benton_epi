# sarsCov2Surveillance.R

# Note: This script does not use the built-in choose(n,k) but instead computes n choose k from the sum/difference of logs.
# This is a bit slower for large n and k, but it can admit much larger n and k then choose(n,k).

# The exact probabiility of observing/encountering x cases given a population = n, sample size k, and count or decimal prevalence prev
sarsCov2Exact <- function(n,    # population
                          k,    # sample size
                          x,    # number of successes
                          prev  # prevalence as decimal or count
                          ) {
  if(prev < 1){prev <- round(n*prev,0)}  # converting decimal prevalence to count
  
  # Handling degenerate cases
  if(min(n,k,x,prev) < 0) {
    print("all arguments must be non-negative")
    return(NULL)
  } else if(k > n) {
    print("k must be less than or equal to n")
    return(NULL)
  } else if(prev > n) {
    print("prev must be less than or equal to n")
    return(NULL)
  } else if(x > n) {
    print("x must be less than or equal to n")
    return(NULL)
  }
  # Handling boundary cases
  else if(k == n){
    if(x == prev) {return(1)}
    else{return(0)}
  } else if(prev == n) {
    if(k == x) {return(1)}
    else{return(0)}
  } else if(k < x | prev < x) {
    return(0)
  } else if(n*k*prev == 0) {
    return(1)
  } else if(k-x > n-prev) {
    return(0)
  }
  # need to break up interior cases because there may be factorial(0) in the formula, which returns an error with log
  else if(x == 0) {
    toplog <- sum(log(c((n-prev-k+1):(n-prev))))-sum(log(c(1:k)))
    botlog <- sum(log(c(n-k+1):n))-sum(log(c(1:k)))
    return(exp(toplog-botlog))
  } else if(k == x) {
    toplog <- sum(log(c((prev-x+1):prev)))-
              sum(log(c(1:x)))
    botlog <- sum(log(c(n-k+1):n))-sum(log(c(1:k)))
    return(exp(toplog-botlog))
  } else {
    toplog <- sum(log(c((prev-x+1):prev)))-
              sum(log(c(1:x)))+
              sum(log(c((n-prev-(k-x)+1):(n-prev))))-
              sum(log(c(1:(k-x))))
    botlog <- sum(log(c(n-k+1):n))-sum(log(c(1:k)))
    return(exp(toplog-botlog))
  }
}

# Probability of success in observing/encountering one or more cases given population = n, sample size = k, prevalence = prev
sarsCov2Success <- function(n,k,prev){
  if(prev < 1){prev <- round(n*prev,0)}  # converting decimal prevalence to count
  
  p <- 0  # initializing probability
  # calls sarsCov2Exact for each case from 1 to total prevalence
  for(x in 1:prev){
    p <- p + sarsCov2Exact(n,k,x,prev)
  }
  return(p)
}

# Generates the mean and standard deviation for the distribution of the
# random variable X = how many times one or more cases observed/encountered out of nSamp,
# with population = n, sample size = k, and prevalence = prev
# The distribution is binomial with parameters nSamp and probability calculated by sarsCov2Success
sarsCov2Distribution <- function(n,    # population
                                 k,    # size of sample
                                 prev, # prevalence
                                 nSamp # number of samples
                                 ){
  p <- sarsCov2Success(n,k,prev)
  sd <- sqrt(nSamp*p*(1-p))
  return(list(mean = nSamp*p,stdev = sd))
}

# Generates random samples drawn from sarsCov2Distribution.
sarsCov2Samples <- function(n,      # population
                            k,      # subpopulation size
                            prev,   # prevalence as count
                            nSamp,  # number of samples in trial
                            nTrials # number of trials
                            ) {
  p <- sarsCov2Success(n,k,prev)
  return(rbinom(nTrials,nSamp,p))
}


# Computes the statistical Power of syndromic surveillance in a subgroup within a community.
# Suppose you know the prevalence in the community and you want to conduct syndromic surveillance
# in a subgroup like a school or long term care facility. What is the likelihood that you will detect
# the disease in the subgroup under the assumption that the disease is present in the subgroup?
# This depends on the likelihood that the disease is in the subgroup, specifically how many cases
# are in the subgroup, which is determined by sarsCov2Exact. It also depends on the likelihood that
# the one or more case will be observed given a certain number of cases within the subgroup. This is
# determined by sarsCov2Success for each possible number of cases in the subgroup.
# The purpose of this function is to give decision makers information about how many people need
# to be tested for adequate syndromic surveillance given a certain community prevalence.
sarsCov2Power <- function(N,              # population of community
                          n,              # population of subgroup
                          k,              # size of sample in subgroup that will be tested
                          prev,           # community prevalence as decimal or count
                          asymp = 1,      # proportion of cases that are asymptomatic
                          sensitivity = 1 # sensitivity of test
                          ) {
  if(prev < 1) {prev <- round(N*prev,0)}           # converts proportion to count
  prev <- ceiling(prev*asymp)                      # surveillance is only for asymptomatic cases
  pwrInverse <- 0                                  # initializing pwrInverse = P(no detection | presence of disease)
  # given prevalence (count) in community, there can be between 0 and min(number of individuals in subgroup,number infections in community) cases.
  # for each possible positive number of cases x, compute P(x|prevalence in community)*P(no detection|x cases in subgroup)
  # sum these probabilities to get pwrInverse
  for(x in 1:min(n,prev)){
    pwrInverse <- pwrInverse + sarsCov2Exact(N,n,x,prev)*(1-sarsCov2Success(n,k,x))
  }
  # power is probability of detecting given presence of disease = (1-P(missing | presence))*(sensitivity of the test)
  pwr <- (1-pwrInverse)*sensitivity                
  return(pwr)
}
