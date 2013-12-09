### StartingPoints ###
# This function will decide the starting points of the 
# Adaptive Rejection Sampling function.
## Inputs: 'D' is an interval specifying the lower and upper bounds of the target density 
##         'h' is the log of the target density
## Output: 'T' is a row vector of two elements


startingpoints <- function(D,h){
#When the user doesn't specify the domain, the input argument could be NA's.
  D[1][which(is.na(D[1])==TRUE)] <- -Inf #turn NA to Infinity
  D[2][which(is.na(D[2])==TRUE)] <- Inf
  if (D[1]==-Inf && D[2]==Inf){
    a <- -4  #we pick -4 and 4 as the default 
    b <- 4
  }
  else if (D[1]==-Inf){
    a <- -4
    b <- D[2]
    if (a >= b){ #check whether the upper bound is lower than -4
      a <- 2*b   #repick a as twice of upper bound
    }
  }
  else if (D[2]==Inf){
    a <- D[1]
    b <- 4
    if (a >= b){ #check whether the lower bound is biger than 4
      b <- 2*a   #repick b as twice of lower bound
    }
  }
  else {
    a <- D[1]
    b <- D[2]
      }
#  return(c(a,b))   
  ap = diag(attr(numericDeriv(quote(h(a)),'a'),'gradient')) #h'(a)
  bp = diag(attr(numericDeriv(quote(h(b)),'b'),'gradient')) #h'(b)
#  return(c(ap,bp)) 
  if (ap > 0 & bp < 0){  #check whether we have an optim between a and b
    T <- c(a,b)   #assign a and b as a row vector
  }
  else {
    stop("Please give a valid domain for g(x)!\nPlease try again!")
  }
  return(T)
}


