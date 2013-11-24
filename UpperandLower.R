### Code to create Upper and Lower Hulls###
createLowHull = function(T, h, D) {
  
  LowerBound = list()
  for (i in 1:(length(T)-1)) {
    LowerBound[i] = list()
    LowerBound[i]$m = (h(T[i+1]) - h(T[i]))/(T[i+1] - T[i])
    LowerBound[i]$b = (T[i+1] * h(T[i]) - T[i] * h(T[i+1]))/(T[i+1] - T[i])
    LowerBound[i]$left = T[i]
    LowerBound[i]$right = T[i+1]
  }
  m = (h(T[-1]) - h(T[-length(T)]))/(T[-1] - T[-length(T)])
  b = (T[-1] * h(T[-length(T)]) - T[-length(T)] * h(T[-1]))/(T[-1] - T[-length(T)])
  left = T[-length(T)]
  right = T[-1]
  LowerBound = data.frame(c(m,b,left,right))
  colnames(LowerBound) = c('m', 'b', 'left', 'right')
  
  return(LowerBound)
}



createUpHull = function(T, h, D) {
  
  
  x = T
  m = numericDeriv(quote(h(x)),'x') # slopes
  b = h(T) - m*T # y-ints
  z = (b[-1] - b[-length(b)])/(m[-length(m)] - m[-1]) # intersection points
  prob0 = exp(b)/m * (exp(m * c(z, D[2])) - exp(m * c(D[1],z))) #possibly unnormalized probabilities
  prob = prob0/sum(prob0) # normalized probabilities
  left = c(D[1],z) # left endpoint for each line segment
  right = c(z, D[2]) # right endpoint for each line segment
  UpBound = data.frame(cbind(m,b,z,prob,left,right))
  colnames(UpBound) = c('m', 'b', 'z', 'prob', 'left','right')
  return(UpBound)
}