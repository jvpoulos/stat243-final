### ARS ###
## The purpose of this function is to perform Adaptive Rejection Sampling.
#  Inputs:
#     1.  'g' is a function handle of the (possibly unnormalized) density that 






ars = function(g, D = c(NA,NA), a=NA, b=NA, n=1) {
  
  h = function(x) log(g(x))
  samp = numeric(n)
  #Check for valid inputs
#   if (is.na(a) & is.na(b)) {
#     T_k = startingpoints(D,h,a,b)
#   } else {
#     stop("Either specify a and b, or leave both as NA")
#   }
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
      Up = createUpHull(T_k,h,D)
      Low = createLowHull(T_k,h,D)
    }
    
    
    
  }
  
  return(samp)
}



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
    a <- D[1]
    b <- 4
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
      a <- D[1] + sqrt(.Machine$double.eps)
      b <- D[2] - sqrt(.Machine$double.eps)
    }
    T <- c(a,b)
  }
  
  return(T)
}



createLowHull = function(T, h, D) {
  
  m = (h(T[-1]) - h(T[-length(T)]))/(T[-1] - T[-length(T)])
  b = (T[-1] * h(T[-length(T)]) - T[-length(T)] * h(T[-1]))/(T[-1] - T[-length(T)])
  left = T[-length(T)]
  right = T[-1]
  LowerBound = data.frame(cbind(m,b,left,right))
  colnames(LowerBound) = c('m', 'b', 'left', 'right')
  
  return(LowerBound)
}


## createUpHull ##
# This function will return the necessary information to find the line segments that
# make up the upper boundary.
## Inputs: 'T' is the set/vector of abscissae
##         'h' is the log of the target density
##         'D' is the vector that specifies the domain
## Output: a dataframe with the following information
##         'm' is the vector of slopes of the line segments between two intersection points
##         'b' is the vector of y-intercepts of the line segments
##         'left' is the vector of left starting points of each line segment
##         'right' is the vector of right ending points of each line segment
##         'prob' is the integral of the exponential function raised to the line segment.
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
      IsInf = TRUE
    } else if (!is.infinite(x) & !is.nan(x) & !is.na(x))
      IsInf = FALSE
    
  }
  
  return(x)
}



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