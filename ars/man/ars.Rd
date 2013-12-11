\name{ars}
\alias{ars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ars(g, D = c(NA, NA), a = NA, b = NA, n = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (g, D = c(NA, NA), a = NA, b = NA, n = 1) 
{
    h = function(x) log(g(x))
    samp = numeric(n)
    T_k = startingpoints(D, h, a, b)
    if (T_k[1] < D[1] | T_k[2] > D[2]) {
        stop("At least one of the starting points is outside the domain.")
    }
    if (g(T_k[1]) <= 0 | g(T_k[2]) <= 0) {
        stop("Invalid starting points")
    }
    Low = createLowHull(T_k, h, D)
    Up = createUpHull(T_k, h, D)
    k = 0
    while (k < n) {
        x.star = sampleUp(Up)
        u = runif(1)
        evals = evalSampPt(x.star, Up, Low)
        testforreject = rejectiontest(x.star, u, evals[1], evals[2], 
            h)
        A = testforreject[1]
        Ups = testforreject[2]
        islogconcave = testforreject[3]
        if (islogconcave == FALSE) {
            stop("Error: Density is not log concave")
        }
        if (A) {
            k = k + 1
            samp[k] = x.star
            if (Ups) {
                T_k = sort(c(T_k, x.star))
            }
        }
        else {
            T_k = sort(c(T_k, x.star))
            Up = createUpHull(T_k, h, D)
            Low = createLowHull(T_k, h, D)
        }
    }
    return(samp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
