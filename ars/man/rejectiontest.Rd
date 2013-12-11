\name{rejectiontest}
\alias{rejectiontest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rejectiontest(x_star, w, l_k, u_k, h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x_star}{
%%     ~~Describe \code{x_star} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}
  \item{l_k}{
%%     ~~Describe \code{l_k} here~~
}
  \item{u_k}{
%%     ~~Describe \code{u_k} here~~
}
  \item{h}{
%%     ~~Describe \code{h} here~~
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
function (x_star, w, l_k, u_k, h) 
{
    if (w <= exp(l_k - u_k)) {
        A = TRUE
        Up = FALSE
        logconcave = TRUE
    }
    else if (w <= exp(h(x_star) - u_k)) {
        A = TRUE
        Up = TRUE
        logconcave = l_k <= u_k
    }
    else {
        A = FALSE
        Up = TRUE
        logconcave = l_k <= u_k
    }
    return(c(A, Up, logconcave))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
