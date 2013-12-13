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

w = runif(0,1)
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

if !logconcave{
  stop('Error: The target density is not log-concave.')
}


