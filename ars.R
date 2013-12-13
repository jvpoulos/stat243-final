#' @name ars
#' 
#' @title Code to implement Adaptive Rejection Sampling (ARS) function
#' 
#' @description
#' This function will perform Adaptive Rejection Sampling
#' 
#' @usage ars(g, D, a = -4, b = 4, n = 1)
#' 
#' @param g The density that we want to sample from (should be a function)
#' @param D A vector specifying the domain on which to sample on.  Defaults to (-Inf, Inf).
#' @param a A possible starting point for the abscissae (will default to -4 if unspecified)
#' @param b A possible starting point for the abscissae (will default to 4 if unspecified)
#' @param n Number of sample points to generate.
#' 
#' @return sample A vector of n elements, each element being a sample from g
#' 
#' @author 
#' Steven Chang \email{stchang22@berkeley.edu}; Willy Lai \email{sv3323@berkeley.edu}; 
#' Jason Poulos \email{poulos@berkeley.edu}; Shengying Wang \email{swang10@berkeley.edu}
#' 
#' @references
#' Gilks, W.R., P. Wild. (1992) Adaptive Rejection Sampling for Gibbs Sampling, Applied Statistics
#' 41:337–348.
#' 
#' @examples
#' g <- function(x) dnorm(x)
#' D <- c(-Inf, Inf)
#' ars(g, D, a = NA, b = NA, n = 10)

ars = function(g, D = c(NA,NA), a=NA, b=NA, n=1) {
  if (is.numeric(D[1]) & is.numeric(D[2])) {
    if (D[1] >= D[2]) {
      stop('Please specify domain in ascending order, i.e. (2,5)')
    }
  }
  h = function(x) log(g(x))
  samp = numeric(n)
  if ((is.na(a) + is.na(b)) == 1) {
    stop('Please either specify both a and b, or do not specify both')
  }
  T_k = startingpoints(D,h,a,b)
  if (T_k[1] < D[1] | T_k[2] > D[2]) {
    stop("At least one of the starting points is outside the domain.")
  }
  if (g(T_k[1]) <= 0 | g(T_k[2]) <=0) {
    stop("Invalid starting points")
  }
  Low = createLowHull(T_k,h,D)
  Up = createUpHull(T_k,h,D)
  k = 0
  while(k < n) {
    x.star = sampleUp(Up)
    u = runif(1)
    evals = evalSampPt(x.star,Up,Low)
    testforreject = rejectiontest(x.star,u,evals[1],evals[2],h)
    A = testforreject[1]
    Ups = testforreject[2]
    islogconcave = testforreject[3]
    if (islogconcave == FALSE) {
      stop("Error: Density is not log concave")
    }
    if (A) {
      k = k + 1
      samp[k] = x.star 
      
      if (Ups) {
        T_k = sort(c(T_k,x.star))
      }
    } else {
      T_k = sort(c(T_k,x.star))
      Up = try(createUpHull(T_k,h,D),silent = TRUE)
      if (is(Up,"try-error")) {
        stop('We may have sampled a number that is too large or too small.\n Please adjust either the starting points or
             the domain.')
      }
      Low = createLowHull(T_k,h,D)
    }
    
  }
  
  return(samp)
}

## startingpoints ##
# This function will decide the starting points of the 
# Adaptive Rejection Sampling function.
## Inputs: 'D' is an interval specifying the lower and upper bounds of the target density 
##         'h' is the log of the target density
## Output: 'T' is a row vector of two elements

startingpoints <- function(D,h,A,B){
  #When the user doesn't specify the domain, the input argument could be NA's.
  D[1][which(is.na(D[1])==TRUE)] <- -Inf #turn NA to Infinity
  D[2][which(is.na(D[2])==TRUE)] <- Inf
  if (D[1]==-Inf && D[2]==Inf){
    if (is.numeric(A) & is.numeric(B)) {
      a <- A
      b <- B
    } else {
      a <- -4  #we pick -4 and 4 as the default 
      b <- 4
    }
    ap = diag(attr(numericDeriv(quote(h(a)),'a'),'gradient')) #h'(a)
    bp = diag(attr(numericDeriv(quote(h(b)),'b'),'gradient')) #h'(b)
    #  return(c(ap,bp)) 
    if (ap > 0 & bp < 0){  #check whether we have an optim between a and b
      T <- c(a,b)   #assign a and b as a row vector
    } else if (ap > 0 & bp >= 0){
      stop("No points to the right of the mode")
    } else if (ap <= 0 & bp < 0) {
      stop("No points to the left of the mode")
    }
    else {
      stop("Please give a valid domain for g(x)!\nPlease try again!")
    }
  }
  else if (D[1]==-Inf){
    a <- -4
    b <- D[2]
    if (a >= b){ #check whether the upper bound is lower than -4
      a <- 2*b   #repick a as twice of upper bound
    }
    T <- c(a,b)
  }
  else if (D[2]==Inf){
    if (is.numeric(A)) {
      a <- A
    } else {
      a <- D[1]
    }
    if (is.numeric(B)) {
      b = B
    } else {
      b = 4
    }
    if (a >= b){ #check whether the lower bound is biger than 4
      b <- 2*a   #repick b as twice of lower bound
    }
    T <- c(a,b)
  }
  else {
    if (is.numeric(A) & is.numeric(B)) {
      a = A
      b = B
    } else {
      a <- D[1] + 0.003
      b <- D[2] - 0.003
    }
    T <- c(a,b)
  }
  
  return(T)
}

#' @include evalSampPt.R
#' @include sampleFromUp.R
#' @include StartingPoints.R
#' @include Rejection_Test.R
#' 
#' @name createLowHull
#' 
#' @title Code to create Upper and Lower Hulls
#' 
#' @description
#' This function will return a dataframe with the information
#' needed to plot the line segments.
#' 
#' @param T set/vector of abscissae
#' @param h log of the target density
#' @param D vector that specifies the domain
#' 
## Output: A dataframe with the following
#' @return m vector of slopes of the line segments between two abscissae
#' @return b vector of y-intercepts of the line segments
#' @return left vector of left starting points of each line segment
#' @return right vector of right ending points of each line segment

createLowHull = function(T, h, D) {
  
  m = (h(T[-1]) - h(T[-length(T)]))/(T[-1] - T[-length(T)])
  b = (T[-1] * h(T[-length(T)]) - T[-length(T)] * h(T[-1]))/(T[-1] - T[-length(T)])
  left = T[-length(T)]
  right = T[-1]
  LowerBound = data.frame(cbind(m,b,left,right))
  colnames(LowerBound) = c('m', 'b', 'left', 'right')
  
  return(LowerBound)
}

#' @name createUpHull 
#' 
#' @description 
#' This function will return the necessary information to find the line segments that
#' make up the upper boundary.
#' 
#' @param T set/vector of abscissae
#' @param h log of the target density
#' @param D vector that specifies the domain
#' 
## Output: a dataframe with the following information
#' @return m vector of slopes of the line segments between two intersection points
#' @return b vector of y-intercepts of the line segments
#' @return left vector of left starting points of each line segment
#' @return right vector of right ending points of each line segment
#' @return prob integral of the exponential function raised to the line segment.

createUpHull = function(T, h, D) {
  
  x = T
  m = diag(attr(numericDeriv(quote(h(x)),'x'),'gradient')) # slopes
  b = h(T) - m*T # y-ints
  z = (b[-1] - b[-length(b)])/(m[-length(m)] - m[-1]) # intersection points
  prob0 = exp(b)/m * (exp(m * c(z, D[2])) - exp(m * c(D[1],z))) #possibly unnormalized probabilities
  prob = prob0/sum(prob0) # normalized probabilities
  prob[is.nan(prob)] = 1
  left = c(D[1],z) # left endpoint for each line segment
  right = c(z, D[2]) # right endpoint for each line segment

  if (m[1] == m[2]) {
    UpBound = data.frame(cbind(m[1], b[1], prob[1], left[1], right[2]))
  } else {
    UpBound = data.frame(cbind(m,b,prob,left,right))
  }
  colnames(UpBound) = c('m', 'b', 'prob', 'left','right')
  return(UpBound)
}


#' @include evalSampPt.R
#' @include UpperandLower.R
#' @include StartingPoints.R
#' @include Rejection_Test.R
#' 
#' @name sampleUp 
#' 
#' @title This function will generate a potential sample point from
#' the enveloping function.
#' 
#' @param UpperHull the dataframe generated by createUpHull
#' 
#' @return x a single point
sampleUp = function(UpperHull) {
  
  emp.cdf = cumsum(UpperHull$prob)
  IsInf = TRUE # repeat until sample is valid
  while(IsInf) {
    u = runif(1)
    ind = min(which(u<emp.cdf, arr.ind = TRUE))
    u = runif(1)
    m = UpperHull$m[ind]
    b = UpperHull$b[ind]
    left = UpperHull$left[ind]
    right = UpperHull$right[ind]
    
    x = log( u*(exp(m*right) - exp(m*left)) + exp(m*left) ) / m
    if (x < UpperHull$left[1] | x > UpperHull$right[length(emp.cdf)]) {
      IsInf = TRUE #make sure the point is within the domain
    } else if (!is.infinite(x) & !is.nan(x))
      IsInf = FALSE
    
  }
  
  return(x)
}

#' @include sampleFromUp.R
#' @include UpperandLower.R
#' @include StartingPoints.R
#' @include Rejection_Test.R
#' 
#' @name evalSampPt
#' 
#' @title Evaluate the enveloping and sandwiching functions
#' 
#' @description
#' This function will evaluate the enveloping and sandwiching functions
#' at the point generated by sampleUp. The returned values are meant
#' to be used to see whether to reject the sample point or not.
#'
#' @param x potential sample point
#' @param UpHull dataframe for the upper boundary
#' @param LowHull dataframe for the lower boundary
#' @return lEval value of the sandwiching function at the sample point
#' @return uEval value of the enveloping function at the sample point


evalSampPt = function(x, UpHull, LowHull) {
  
  
  # Evaluate lower bound at sample point
  if (x < min(LowHull$left)) {
    lEval = -Inf
  } else if (x > max(LowHull$right)) {
    lEval = -Inf
  } else {
    # Determine which set of information to use
    ind = which(x >= LowHull$left & x <= LowHull$right, arr.ind = TRUE) 
    lEval = LowHull$m[ind] * x + LowHull$b[ind]
  }
  
  # Evaluate upper bound at sample point
  indR = which(x >= UpHull$left & x <= UpHull$right, arr.ind = TRUE)
  uEval = UpHull$m[indR] * x + UpHull$b[indR]
  
  return(c(lEval, uEval))
}

#' @include evalSampPt.R
#' @include sampleFromUp.R
#' @include UpperandLower.R
#' @include StartingPoints.R
#'
#' @name rejectiontest
#' 
#' @title This code will perform the rejection algorithm.
#' 
#' @description
#'The function rejectiontest returns 3 logicals of whether we accept the sample,
#'have to update h(x), and checks if the function h(x) is log-concave at x_star.  
#'If it is not log-concave, then the function is supposed to stop.
#'
#' @param x_star sample candidate point
#' @param w a number from the Uniform(0,1) distribution
#' @param l_k evaluation of lower bound at x_star
#' @param u_k evaluation of upper bound at x_star
#' @param h log of target density
#' 
#' Output is a logical vector of three elements
#' @return A TRUE if we accept x_star
#' @return Up TRUE if we need to update the abscissae
#' @return logconcave TRUE if log-concave at x_star

rejectiontest = function(x_star,w, l_k, u_k, h){
  if( w <= exp(l_k-u_k)){
    A = TRUE
    Up = FALSE
    logconcave = TRUE
  }
  else if (w <= exp(h(x_star)-u_k)){ 
    A = TRUE
    Up = TRUE
    logconcave = l_k <= u_k
  }
  else{
    A = FALSE
    Up = TRUE
    logconcave = l_k <= u_k
  }
  return(c(A,Up,logconcave))
}


#test for the mainfunction

#test the main function
# g is the sampling density
# D is the domain
# n is the number of points
# p is the pdf of the true sampling density or generated data from a truncated distribution


# We are using the kolmogorov-smirnov test for comparison of sampled points with true distribution
# If the p-value is bigger than 5%, we pass.

#define the testing function
test = function() {
  set.seed(0)
  par(mfrow=c(2,3))
  titlebank = c("Standard Normal", 'Bounded Beta(2,5)', 'Exponential(1.5)', 'Bounded Gamma(2,2)', 'Bounded Weibull(1.5,1)', 'Bounded Logis(5,2)')
  testMain <- function(g,D,a,b,n,p,i){
    x <- ars(g, D, a, b, n)
    hist(x, freq = FALSE, main = paste(c(titlebank[i], "Distribution", sep = " ")))
    curve(g, add = TRUE, col = "blue")
    pvalue <- ks.test(x,p)$p.value
    if (pvalue > 0.05){
      print("Pass the Kolmogorov-Smirnov test at level = 0.05.")
    }
    else { print ("Fail the Kolmogorov-Smirnov test at level = 0.05.")}
  }
  
  #test for standard normal density
  print('Test a standard normal.')
  i = 1
  testMain(dnorm, D = c(-Inf,Inf), a=NA, b=NA, n=10000, pnorm,i)
  
  
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
  print('Test for Bounded Beta(2,5)')
  i = i+1
  testMain(dbetatest, D = c(0.1,0.9), a = 0.15, b = 0.85,n=10000, x,i)
  
  
  #test for exponential with lambda=1.5
  #define the exponential density with lambda=1.5
  dexptest <- function(x){
    y <- dexp(x,1.5)
    return(y)
  }
  pexptest <- function(x) {
    y <- pexp(x,1.5)
    return(y)
  }
  
  print('Test an Exponential with Rate 1.5')
  i = i+1
  testMain(dexptest, D = c(0,Inf), a=NA, b=NA, n=10000, pexptest,i)
  
  
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
  i = i+1
  testMain(dgammatest, D = c(0,10), a=NA, b=NA, n=10000, x,i)
  
  
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
  i = i+1
  testMain(dweibulltest, D = c(0,15), a=NA, b=NA, n=10000, x, i)
  
  
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
  print('Test for Bounded Logis(5,2)')
  i = i+1
  testMain(dlogistest, D = c(-5,15), a=NA, b=NA, n=10000, x,i)
  
}

