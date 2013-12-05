#' @include evalSampPt.R
#' @include sampleFromUp.R
#' 
#' @name createLowHull
#' 
#' @title Code to create Upper and Lower Hulls###
#' @description
#' This function will return a dataframe with the information
#' needed to plot the line segments.
#' @param T set/vector of abscissae
#' @param h log of the target density
#' @param D vector that specifies the domain
## Output: A dataframe with the following
#' @return m is the vector of slopes of the line segments between two abscissae
#' @return b is the vector of y-intercepts of the line segments
#' @return left is the vector of left starting points of each line segment
#' @return right is the vector of right ending points of each line segment

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
  UpBound = data.frame(cbind(m,b,prob,left,right))
  colnames(UpBound) = c('m', 'b', 'prob', 'left','right')
  return(UpBound)
}