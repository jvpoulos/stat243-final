sampleUp = function(UpperHull) {
  
  emp.cdf = cumsum(UpperHull$prob)
  IsInf = TRUE # repeat until sample is valid
  while(IsInf) {
    u = rnorm(1)
    ind = min(which(u<emp.cdf, arr.ind = TRUE))
    u = rnorm(1)
    m = UpperHull$m[ind]
    b = UpperHull$b[ind]
    left = UpperHull$left[ind]
    right = UpperHull$right[ind]
  
    x = log( U*(exp(m*right) - exp(m*left)) + exp(m*left) ) / m
    if (!is.infinite(x) & !is.nan(x)) {
      IsInf = FALSE
    }
  
  
  }
  
  return(x)
}