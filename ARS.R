#' A package to implement adaptive-rejection sampler
#' @name ARS
#' @docType package
#' 
#' 

ARS = function(g, D, a=NA, b=NA, n=1) {
  
  h = function(x) log(g(x))
  samp = numeric(n)
  #Check for valid inputs
  if (is.na(a) & is.na(b)) {
    T_k = StartingPoints(D,h)
  } else if (sum(is.na(c(a,b))) == 0 & diag(attr(numericDeriv(quote(h(a)),'a'))) < 0 & diag(attr(numericDeriv(quote(h(b)),'b') > 0))) {
    T_k = c(a,b)
  }  else {
    stop("Either specify a and b, or leave both as NA")
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
    if (islogconcave == TRUE) {
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
  
  
} 