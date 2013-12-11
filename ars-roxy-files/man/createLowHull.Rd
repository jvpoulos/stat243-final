\name{createLowHull}
\alias{createLowHull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
createLowHull(T, h, D)
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
    m = (h(T[-1]) - h(T[-length(T)]))/(T[-1] - T[-length(T)])
    b = (T[-1] * h(T[-length(T)]) - T[-length(T)] * h(T[-1]))/(T[-1] - 
        T[-length(T)])
    left = T[-length(T)]
    right = T[-1]
    LowerBound = data.frame(cbind(m, b, left, right))
    colnames(LowerBound) = c("m", "b", "left", "right")
    return(LowerBound)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
