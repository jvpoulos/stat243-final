#test for starting points

#test the startingpoints function
# D is the domain
# h is the log of sampling density
# T is the correct starting points
# We are going to compare our starting points from startingpoints function and the correct ones
# If they are the same, we pass.

#define the testing function
testStartingpoints <- function(D,h,A,B,T){
  Tp <- startingpoints(D,h,A,B)
  if (Tp[1] == T[1] && Tp[2] == T[2]){
    print('Pass! The output values are correct.')
  }
  else print('Fail! The output values are incorrect.')
}

#define log normal
lognor <- function(x){
  y <- dnorm(x,log=TRUE)
}
#define log beta
logbet <- function(x){
  y <- dbeta(x,2,3,log=TRUE)
  return(y)
}
#define log exp
logexp <- function(x){
  y <- dbeta(x,log=TRUE)
  return(y)
}

#test for multiple cases
testStartingpoints(D=c(-Inf,Inf),lognor,A=NA,B=NA,T=c(-4,4))
testStartingpoints(D=c(-Inf,Inf),lognor,A=-6,B=5,T=c(-6,5))
testStartingpoints(D=c(0.1,Inf),logexp,A=NA,B=NA,T=c(0.1,4))
testStartingpoints(D=c(-Inf,0.9),lognor,A=NA,B=NA,T=c(-4,0.9))
testStartingpoints(D=c(0.1,0.9),logbet,A=NA,B=NA,T=c(0.1+0.003,0.9-0.003))
testStartingpoints(D=c(0.1,0.9),logbet,A=0.2,B=0.8,T=c(0.2,0.8))