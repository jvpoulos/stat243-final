#This is a rough draft of the code for the rejection sampling test.

#The function rejectiontest returns 3 logicals of whether we accept the sample,
#have to update h(x), and checks if the function h(x) is log-concave at x_star.  
#If it is not log-concave, then the function is supposed to stop.

#Suppose x_star is a sample from s_k

w = runif(0,1)
rejectiontest = function(x_star,w){
  if( w <= exp(l_k(x_star)-u_k(x_star))){
    A = TRUE
    Up = FALSE
    logconcave = TRUE
  }
  elseif (w <= exp(h(x_star)-u_k(x_star)){ 
    A = TRUE
    Up = TRUE
    logconcave = l_k(x_star) <= h(x_star) & h(x_star) <= u_k(x_star)
  }
  else{
    A = FALSE
    Up = TRUE
    logconcave = l_k(x_star) <= h(x_star) & h(x_star) <= u_k(x_star)
  }
  return(c(A,Up,logconcave))
}

if !logconcave{
  stop('Error: The target density is not log-concave.')
}


