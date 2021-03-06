# syndSurveillance.R

# Note: This script does not use the built-in choose(n,k) but instead computes n choose k from the sum/difference of logs.
# This is a bit slower for large n and k, but it can admit much larger n and k then choose(n,k).

# The exact probability of detecting or being exposed to x cases given
# a population = n, sample size k, and count or decimal prevalence prev, with test sensitivity/transmissibilty = sensitivity.
syndExact <- function(n,    # population
                      k,    # sample size
                      x,    # number of successes
                      prev,  # prevalence as decimal or count
                      sensitivity = 1 # sensitivity of the test/likelihood of infection from an encounter
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
  else if(k < x | prev < x) {
    return(0)
  }
  else if(n*k*prev == 0) {
    return(1)
  }
  # need to break up interior cases because there may be factorial(0) in the formula, which returns an error with log
  else if(x == 0) {
    botlog <- sum(log(c(n-k+1):n))-sum(log(c(1:k)))
    p <- 0
    # handle j = 0
    if(k <= n - prev) {
      toplog <- sum(log((n-prev-k+1):(n-prev)))-
        sum(log(1:k))
      p <- p + exp(toplog-botlog)
    }
    for(j in 1:min(k,prev)){
      if(k-j == n-prev | k-j == 0) {
        toplog <- sum(log((prev-j+1):prev))-
          sum(log(1:j))+
          log(dbinom(x,j,sensitivity))
        p <- p + exp(toplog-botlog)
      } else if(k-j < n-prev) {
        toplog <- sum(log((prev-j+1):prev))-
          sum(log(1:j))+
          sum(log((n-prev-(k-j)+1):(n-prev)))-
          sum(log(1:(k-j)))+
          log(dbinom(x,j,sensitivity))
        p <- p + exp(toplog-botlog)
      }
    }
    return(p)
  } else if(k == x) {
    botlog <- sum(log(c(n-k+1):n))-sum(log(c(1:k)))
    toplog <- sum(log(c((prev-x+1):prev)))-sum(log(c(1:x)))+log(sensitivity^x)
    return(exp(toplog-botlog))
  } else {
    botlog <- sum(log(c(n-k+1):n))-sum(log(c(1:k)))
    p <- 0
    for(j in x:min(k,prev)) {
      if(k-j == n-prev | k-j == 0) {
        toplog <- sum(log((prev-j+1):prev))-
          sum(log(1:j))+
          log(dbinom(x,j,sensitivity))
        p <- p + exp(toplog-botlog)
      } else if(k-j < n-prev) {
        toplog <- sum(log((prev-j+1):prev))-
          sum(log(1:j))+
          sum(log((n-prev-(k-j)+1):(n-prev)))-
          sum(log(1:(k-j)))+
          log(dbinom(x,j,sensitivity))
        p <- p + exp(toplog-botlog)
      }
    }
    return(p)
  }
}

# Probability of success in observing/encountering one or more cases given population = n, sample size = k, prevalence = prev
syndSuccess <- function(n, # population
                        k, # sample size
                        prev, # prevalence in population
                        sensitivity=1 # sensitivity of the test/likelihood of infection
){
  if(prev < 1){prev <- round(n*prev,0)}  # converting decimal prevalence to count
  
  p <- 0  # initializing probability
  # calls syndExact for each case from 1 to sample size or total prevalence
  for(x in 1:min(k,prev)){
    p <- p + syndExact(n,k,x,prev,sensitivity)
  }
  return(p)
}

# Generates the mean and standard deviation for the distribution of the
# random variable X = how many times one or more cases observed/encountered out of nSamp,
# with population = n, sample size = k, and prevalence = prev
# The distribution is binomial with parameters nSamp and probability calculated by syndSuccess
syndDistribution <- function(n,    # population
                             k,    # size of sample
                             prev, # prevalence
                             nSamp, # number of samples
                             sensitivity=1 # sensitivity of the test/likelihood of infection
){
  p <- syndSuccess(n,k,prev,sensitivity)
  sd <- sqrt(nSamp*p*(1-p))
  return(list(mean = nSamp*p,stdev = sd))
}

# Generates random samples drawn from syndDistribution.
syndSamples <- function(n,      # population
                        k,      # subpopulation size
                        prev,   # prevalence as count
                        nSamp,  # number of samples in trial
                        nTrials, # number of trials
                        sensitivity=1 # sensitivity of the test/likelihood of infection
) {
  p <- syndSuccess(n,k,prev,sensitivity)
  return(rbinom(nTrials,nSamp,p))
}

# Based on estimated community prevalence, what is the expected number of cases
# that will be identified when testing a sample of size k
syndExpected <- function(n,k,prev,sensitivity) {
  if(prev < 1) {prev <- round(n*prev,0)}           # converts proportion to count
  
  expected <- 0 # initialize expected value
  for(x in 1:min(k,prev)) {
    expected <- expected +
      x * syndExact(n,k,x,prev,sensitivity)
  }
  return(expected)
}

# Computes something analogous to statistical Power of syndromic surveillance in a subgroup within a community.
# Suppose you know the prevalence in the community and you want to conduct syndromic surveillance
# in a subgroup like a school or long term care facility. What is the likelihood that you will detect
# the disease in the subgroup given the prevalence in the community? Or more importantly, what is the likelihood
# that you will not miss the presence of the disease?
# This depends on the likelihood that the disease is in the subgroup, specifically how many cases
# are in the subgroup, which is determined by syndExact. It also depends on the likelihood that
# the one or more case will be observed given a certain number of cases within the subgroup. This is
# determined by syndSuccess for each possible number of cases in the subgroup.
# The purpose of this function is to give decision makers information about how many people need
# to be tested for adequate syndromic surveillance given a certain community prevalence.
# If k=0 people are tested, things get a little funky, since there is still the probability that
# no cases are in the subgroup n. So syndPower(N,n,0,prev) = P(0 cases in subgroup n | prevalence in N)
# instead of the Power of detection.
syndPower <- function(N,              # population of community
                      n,              # population of subgroup
                      k,              # size of sample in subgroup that will be tested
                      prev,           # community prevalence as decimal or count
                      infLikelihood = 1, # Likelihood of infection
                      sensitivity = 1, # sensitivity of test
                      asymp = 1      # proportion of cases that are asymptomatic
) {
  if(prev < 1) {prev <- round(N*prev,0)}           # converts proportion to count
  prev <- ceiling(prev*asymp)                      # surveillance is only for asymptomatic cases
  pwrInverse <- 0                                  # initializing pwrInverse = P(no detection | presence of disease)
  # given prevalence (count) in community, there can be between 0 and min(number of individuals in subgroup,number infections in community) cases.
  # for each possible positive number of cases x, compute P(x|prevalence in community)*P(no detection|x cases in subgroup)
  # sum these probabilities to get pwrInverse
  for(x in 1:min(n,prev)){
    pwrInverse <- pwrInverse + syndExact(N,n,x,prev,sensitivity = infLikelihood)*(1-syndSuccess(n,k,x,sensitivity = 1))
  }
  # power is probability of detecting given presence of disease = (1-P(missing | presence))
  pwr <- (1-pwrInverse)                
  return(pwr)
}
