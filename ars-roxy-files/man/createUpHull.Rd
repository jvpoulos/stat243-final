\name{createUpHull}
\alias{createUpHull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
createUpHull(T, h, D)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{
%%     ~~Describe \code{T} here~~
}
  \item{h}{
%%     ~~Describe \code{h} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (T, h, D) 
{
    x = T
    m = diag(attr(numericDeriv(quote(h(x)), "x"), "gradient"))
    b = h(T) - m * T
    z = (b[-1] - b[-length(b)])/(m[-length(m)] - m[-1])
    prob0 = exp(b)/m * (exp(m * c(z, D[2])) - exp(m * c(D[1], 
        z)))
    prob = prob0/sum(prob0)
    prob[is.nan(prob)] = 1
    left = c(D[1], z)
    right = c(z, D[2])
    UpBound = data.frame(cbind(m, b, prob, left, right))
    colnames(UpBound) = c("m", "b", "prob", "left", "right")
    return(UpBound)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
