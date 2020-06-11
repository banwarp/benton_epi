### Readme file for syndSurveillance
This file contains descriptions of the functions and an example.

#### Please check my work - I think I have the counting correct, but I would appreciate confirmation/corrections.

### Scripts included in [syndSurveillance.R](sundSurveillance.R)
- `syndExact(n,k,x,prev,sensitivity)`
- `syndSuccess(n,k,prev,sensitivity)`
- `syndDistribution(n,k,prev,nSamp,sensitivity)`
- `syndSamples(n,k,prev,nSamp,nTrials,sensitivity)`
- `syndPower(N,n,k,prev,sensitivity,sensitivity)`

### syndExact
This is the base function in syndSurveillance. All other functions use it in some way.  
`syndExact(n,k,x,prev,sensitivity)` computes the probability of observing `x` cases of disease, given a population `n`, a sample size `k`, and a population prevalence (count, not proportion) `prev`. `syndExact` uses a combinatorial counting method to calculate exact probabilities. However, because populations and sample sizes can be very large, using the built-in `choose(n,k)` function can lead to `Inf`. Therefore `syndExact` uses a "log and sum" approach to the combinatorial counting instead of the factorial approach normally used in combinatorics.  Furthermore, when `syndExact` is used for testing or for the probability of infection, `sensitivity` can be used to represent the sensitivity of the test or the likelihood of infection based on contact.

##### Combinatorial method used in `syndExact`
Given a population `n` and a sample size `k`, there are `choose(n,k)` possible samples. Of these samples, we want just those that include exactly `x` cases of disease. This number depends on the prevalence `prev` in the population. Specifically:  

There are `choose(prev,x)` ways to select exactly `x` cases from the total cases in the population. These fill out the first `x` spots in the sample of size `k`. Then there are `choose(n-prev,k-x)` ways to select the remaining `k-x` spots from among the non-infected in the population. Therefore, the number of samples of size `k` with exactly `x` cases is given by:
```
choose(prev,x) * choose(n-prev,k-x)
```
and the probability of exactly `x` cases is given by:
```
p <- choose(prev,x) * choose(n-prev,k-x) / choose(n,k)
```

This the basic formula. However, with large `n`, the factorials involved tend toward infinity. To get around this, I take the logs of the combinatorial formula and get:
log(choose(m,n)) = sum_\[i from 1 to m \] log(i) - sum_\[j from 1 to n\] Log(j) - sum_\[k from 1 to (m-n)\] Log(k)

which simplifies to:

sum_\[i from n+1 to m\] Log(i) - sum_\[j from 1 to n\] Log(j)  

Therefore, in the R code the full equation, along with sensitivity, is:
```
p <- exp(sum(log(c((prev-x+1):prev)))-
         sum(log(c(1:x)))+
         sum(log(c((n-prev-(k-x)+1):(n-prev))))-
         sum(log(c(1:(k-x))))-
         (sum(log(c(n-k+1):n))-
         sum(log(c(1:k))))) * sensitivity^x
```
However, there are degenerate and boundary cases that must be handled. For example, `k > n` is a degenerate case and `k = 0` is a boundary case. Therefore the first half of the function is just for handling those cases:
```
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
  else if(k == n){
    if(x == prev) {return(sensitivity^x)}
    else{return(0)}
  } else if(prev == n) {
    if(k == x) {return(sensitivity^x)}
    else{return(0)}
  } else if(k < x | prev < x) {
    return(0)
  } else if(n*k*prev == 0) {
    return(1)
  } else if(k-x > n-prev) {
    return(0)
  }
```

Furthermore, even the interior cases could result in the `log(0)` appearing in the function, so those cases are broken out as well:
```
# need to break up interior cases because there may be factorial(0) in the formula, which returns an error with log
  else if(x == 0) {
    toplog <- sum(log(c((n-prev-k+1):(n-prev))))-sum(log(c(1:k)))
    botlog <- sum(log(c(n-k+1):n))-sum(log(c(1:k)))
    return(exp(toplog-botlog)*sensitivity^x)
  } else if(k == x) {
    toplog <- sum(log(c((prev-x+1):prev)))-
              sum(log(c(1:x)))
    botlog <- sum(log(c(n-k+1):n))-sum(log(c(1:k)))
    return(exp(toplog-botlog)*sensitivity^x)
  } else {
    toplog <- sum(log(c((prev-x+1):prev)))-
              sum(log(c(1:x)))+
              sum(log(c((n-prev-(k-x)+1):(n-prev))))-
              sum(log(c(1:(k-x))))
    botlog <- sum(log(c(n-k+1):n))-sum(log(c(1:k)))
    return(exp(toplog-botlog)*sensitivity^x)
  }
}
```
The `sensitivity` parameter allows syndExact to adjust to probability for the sensitivity of the test (if using syndExact for testing), or the likelihood of infection (if using syndExact for transmission probabilities).


Here are two examples of comparing the built in `choose(m,n)` and `syndExact`:
```
# Example 1:
n <- 100
k <- 10
x <- 2
prev <- 4
choose(prev,x)*choose(n-prev,k-x)/choose(n,k)
syndExact(n,k,x,prev,sensitivity = 1)

# Example 2:
n <- 400000
k <- 1000
x <- 10
prev <- 800
choose(prev,x)*choose(n-prev,k-x)/choose(n,k)
syndExact(n,k,x,prev,sensitivity = 1)
```

### syndSuccess
The main goal of syndromic surveillance in a setting like the COVID-19 pandemic is to determine if the disease exists in the sample at all. It is not as important exactly how many cases, just if there is a positive number. `syndSuccess(n,k,prev)` computes the probability of observing at least 1 case given a population `p`, sample size `k`, and prevalence `prev`. It works by calling `syndExact` for all `x` from `1` to `prev` (or `prev*n` if `prev` is decimal). `syndSuccess` passes the sensitivity of the test/likelihood of transmission to `syndExact`.
```
syndSuccess <- function(n, # population
                        k, # sample size
                        prev, # prevalence in population
                        sensitivity=1 # sensitivity of the test/likelihood of infection
                        ){
  if(prev < 1){prev <- round(n*prev,0)}  # converting decimal prevalence to count
  
  p <- 0  # initializing probability
  # calls syndExact for each case from 1 to total prevalence
  for(x in 1:prev){
    p <- p + syndExact(n,k,x,prev,sensitivity)
  }
  return(p)
}
```

### syndDistribution and syndSamples
Although syndromic surveillance is not as concerned with a distribution or generating samples, other thought experiments do rely on these functions. For example, if you meet 50 random people out of a population of 1,000 every day for 30 days, and the population prevalence is 2 cases, on how many days would you expect to have at least one exposure? This is answered with `syndDistribution(1000,50,2,30)`, which returns the expected value (mean) and the standard deviation. `syndDistribution` calls `syndSuccess` and then uses the generated probability in a binomial distribution with size `nSamp`. `syndDistribution` passes the sensitivity of the test/likelihood of transmission to `syndExact`.
```
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
```
`syndSamples` takes this one step further and samples `syndDistribution` to generate a random variable from the distribution.
```
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
```

### syndPower
Along with `syndExact`, `syndPower` is the function that requires the closest inspection. The goal is to compute an analogue of statistical Power for syndromic surveillance in a subgroup within a community where the disease exists. Suppose you know the prevalence `prev` in the community with population `N` and you want to conduct syndromic surveillance in a subgroup like a school or long term care facility with population `n`. What is the likelihood that you will detect the disease in the subgroup by sampling `k` people? This depends on the likelihood that the disease is in the subgroup, specifically how many cases are in the subgroup, which is determined by `syndExact`. It also depends on the likelihood that the one or more case will be observed given a certain number of cases within the subgroup. This is determined by `syndSuccess` for each possible number of cases in the subgroup. In order to compute the Power, we sum the product of `syndExact(N,n,prev) * (1-synSuccess(n,k,x))` for every `x` from 1 to the minimum of `prev` and `k`. This gives us the probability of failing to detect the disease given its presence, so Power is the additive complement (i.e. 1-P(failing to detect | presence of disease)).  
Furthermore, detection of disease depends on the sensitivity of the test, so `sensitivity` is passed to `syndExact` and `syndSuccess`:
Also, if we assume that any symptomatic person will be tested immediately, then we really only care about syndromic surveillance if the whole subgroup is asymptomatic. `asymp` reduces the prevalence to only asymptomatic people. E.g. `asymp = .35` means that 35% of cases are asymptomatic, and `syndPower` uses `prev*asymp` as its prevalence.  

```
syndPower <- function(N,              # population of community
                          n,              # population of subgroup
                          k,              # size of sample in subgroup that will be tested
                          prev,           # community prevalence as decimal or count
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
    pwrInverse <- pwrInverse + syndExact(N,n,x,prev,sensitivity)*(1-syndSuccess(n,k,x,sensitivity))
  }
  # power is probability of detecting given presence of disease = (1-P(missing | presence))
  pwr <- (1-pwrInverse)                
  return(pwr)
}
```

#### Example of using `syndPower`
Suppose that the community population is `N=8`, the subgroup population is `n=5`, you will sample `k=3` people, and the community prevalence of asymptomatic cases is `prev=2`.  
- There are `choose(N,n) = choose(8,5)` possible subgroups.
```
combn(8,5)
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18]
[1,]    1    1    1    1    1    1    1    1    1     1     1     1     1     1     1     1     1     1
[2,]    2    2    2    2    2    2    2    2    2     2     2     2     2     2     2     2     2     2
[3,]    3    3    3    3    3    3    3    3    3     3     4     4     4     4     4     4     5     5
[4,]    4    4    4    4    5    5    5    6    6     7     5     5     5     6     6     7     6     6
[5,]    5    6    7    8    6    7    8    7    8     8     6     7     8     7     8     8     7     8
     [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35]
[1,]     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
[2,]     2     2     3     3     3     3     3     3     3     3     3     3     4     4     4     4     5
[3,]     5     6     4     4     4     4     4     4     5     5     5     6     5     5     5     6     6
[4,]     7     7     5     5     5     6     6     7     6     6     7     7     6     6     7     7     7
[5,]     8     8     6     7     8     7     8     8     7     8     8     8     7     8     8     8     8
     [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46] [,47] [,48] [,49] [,50] [,51] [,52]
[1,]     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     3     3
[2,]     3     3     3     3     3     3     3     3     3     3     4     4     4     4     5     4     4
[3,]     4     4     4     4     4     4     5     5     5     6     5     5     5     6     6     5     5
[4,]     5     5     5     6     6     7     6     6     7     7     6     6     7     7     7     6     6
[5,]     6     7     8     7     8     8     7     8     8     8     7     8     8     8     8     7     8
     [,53] [,54] [,55] [,56]
[1,]     3     3     3     4
[2,]     4     4     5     5
[3,]     5     6     6     6
[4,]     7     7     7     7
[5,]     8     8     8     8
```
#### Selecting the subgroup
From a combinatorial perspective, it doesn't matter whether the subgroup individuals have fixed numbers (e.g. 1,2,3,4,5) and the cases are randomly assigned to them, or whether the cases have fixed numbers (e.g. 1,2) and they are randomly assigned to the subgroup. We use the second approach. So assume that individuals 1 and 2 are the cases. 30 of the 56 combinations include exactly 1 case, and a separate 20 include exactly two cases. These 50 combinations represent the state of the world where there is one or more diseases in the subgroup.  

We want to know the probability that we fail to observe the disease if there is exactly 1 case, and if there are exactly 2 cases in our subgroup.  

##### Exactly 1 case in our subgroup
There are `choose(n,k) = choose(5,3)` ways to select our sample for testing. By renumbering, assume individual 1 is the case.
```
combn(5,3)
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
[1,]    1    1    1    1    1    1    2    2    2     3
[2,]    2    2    2    3    3    4    3    3    4     4
[3,]    3    4    5    4    5    5    4    5    5     5
```
4 of the 10 combinations exclude individual 1.

##### Exactly 2 cases in our subgroup
There are still `choose(5,3)` ways to select our sample for testing, and now individuals 1 and 2 are cases.  Now only 1 of the 10 combinations excludes both individuals 1 and 2.

##### Computing Power
The probability that we fail to detect the disease given that it is present is therefore ` (30/56) * (4/10) + (20/56) * (1/10) = 0.25`. If we assume the sensitivity of the test is 100%, the Power of our sampling scheme is `(1-0.25)*1 = 0.75`. You can confirm that `syndPower(8,5,3,2,sensitivity=1,asymp=1)` = 0.75
