\name{sampleUp}
\alias{sampleUp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleUp(UpperHull)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UpperHull}{
%%     ~~Describe \code{UpperHull} here~~
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
function (UpperHull) 
{
    emp.cdf = cumsum(UpperHull$prob)
    IsInf = TRUE
    while (IsInf) {
        u = runif(1)
        ind = min(which(u < emp.cdf, arr.ind = TRUE))
        u = runif(1)
        m = UpperHull$m[ind]
        b = UpperHull$b[ind]
        left = UpperHull$left[ind]
        right = UpperHull$right[ind]
        x = log(u * (exp(m * right) - exp(m * left)) + exp(m * 
            left))/m
        if (!is.infinite(x) & !is.nan(x) & !is.na(x)) {
            IsInf = FALSE
        }
    }
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
