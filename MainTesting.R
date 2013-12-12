#test for the mainfunction

#test the main function
# g is the sampling density
# D is the domain
# n is the number of points
# p is the pdf of the true sampling density or generated data from a truncated distribution


# We are using the kolmogorov-smirnov test for comparison of sampled points with true distribution
# If the p-value is bigger than 5%, we pass.

#define the testing function
BigTestFunction = function() {
  set.seed(0)
  testMain <- function(g,D,a,b,n,p){
    x <- ars(g, D, a, b, n)
    hist(x, freq = FALSE)
    curve(g, add = TRUE, col = "blue")
    pvalue <- ks.test(x,p)$p.value
    if (pvalue > 0.05){
      print("Pass the Kolmogorov-Smirnov test at level = 0.05.")
    }
    else { print ("Fail the Kolmogorov-Smirnov test at level = 0.05.")}
  }

  #test for standard normal density
  print('Test a standard normal.')
  testMain(dnorm, D = c(-Inf,Inf), a=NA, b=NA, n=10000, pnorm)

  #test for truncated beta(2,5)
  #define beta(2,5) density
  dbetatest <- function(x){
    y <- dbeta(x,2,5)
    return(y)
  }

  #function for sampling from truncated beta distribution
  rtrunbeta <- function(n, shape1, shape2, lower, upper){
    y <- rep(NA, length = n)
    k <- 0
    while(k < n){
      cand <- rbeta(1, shape1, shape2)
      if(cand > lower & cand < upper){
        k <- k+1
        y[k] <- cand
      }
    }
    return(y)
  }
  #generate 10000 random truncated beta(2,5) in [0.1,0.9]
  x <- rtrunbeta(10000,2,5,0.1,0.9) 
  #test for truncated beta(2,5) in [0.1,0.9]
  print('Test for Beta(2,5)')
  testMain(dbetatest, D = c(0.1,0.9), a = 0.15, b = 0.85,n=10000, x)

  #test for exponential with lambda=1.5
  #define the exponential density with lambda=1.5
  dexptest <- function(x){
    y <- dexp(x,1.5)
    return(y)
  }
  #define the exponential pdf with lambda=1.5
  pexptest <- function(x){
    y <- pexp(x,1.5)
    return(y)
  }
  print('Test an Exponential with Rate 1.5')
  testMain(dexptest, D = c(0,Inf), a=NA, b=NA, n=10000, pexptest)
  
  #test for truncated gamma(2,2)
  #define gamma(2,2) density
  dgammatest <- function(x){
    y <- dgamma(x,shape=2, rate=2)
    return(y)
  }
  #function for sampling from truncated gamma distribution
  rtrungamma <- function(n, shape, rate, lower, upper){
    y <- rep(NA, length = n)
    k <- 0
    while(k < n){
      cand <- rgamma(1, shape, rate)
      if(cand > lower & cand < upper){
        k <- k+1
        y[k] <- cand
      }
    }
    return(y)
  }
  #generate 10000 random truncated gamma(2,2) in [0,10]
  x <- rtrungamma(10000,2,2,0,10) 
  #test for truncated gamma(2,2) in [0,10]
  print('Test for Bounded Gamma(2,2)')
  testMain(dgammatest, D = c(0,10), a=NA, b=NA, n=10000, x)
  
  #test for truncated weibull density
  #define weibull(1.5,1) density
  dweibulltest <- function(x){
    y <- dweibull(x,shape=1.5, scale=1)
    return(y)
  }

  #function for sampling from truncated weibull distribution
  rtrunweibull <- function(n, shape, scale, lower, upper){
    y <- rep(NA, length = n)
    k <- 0
    while(k < n){
      cand <- rweibull(1, shape, scale)
      if(cand > lower & cand < upper){
        k <- k+1
        y[k] <- cand
      }
    }
    return(y)
  }
  #generate 10000 random truncated weibull(1.5,1) in [0,15]
  x <- rtrunweibull(10000,1.5,1,0,15) 
  #test for truncated weibull(1.5,1) in [0,15]
  print('Test for Bounded Weibull(1.5,1)')
  testMain(dweibulltest, D = c(0,15), a=NA, b=NA, n=10000, x)

  #test for truncated logistic density
  #define logis(5,2) density
  dlogistest <- function(x){
    y <- dlogis(x,location = 5, scale=2)
    return(y)
  }

  #function for sampling from truncated logistic distribution
  rtrunlogis <- function(n, location, scale, lower, upper){
    y <- rep(NA, length = n)
    k <- 0
    while(k < n){
      cand <- rlogis(1, location, scale)
      if(cand > lower & cand < upper){
        k <- k+1
        y[k] <- cand
      }
    }
    return(y)
  }
  #generate 10000 random truncated logis(5,2) in [-5,15]
  x <- rtrunlogis(10000,5,2,-5,15) 
  #test for truncated logis(5,2) in [-5,15]
  print('Test for bounded Logis(5,2)')
  testMain(dlogistest, D = c(-5,15), a=NA, b=NA, n=10000, x)
}
