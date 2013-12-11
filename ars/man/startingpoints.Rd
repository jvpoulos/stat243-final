\name{startingpoints}
\alias{startingpoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
startingpoints(D, h, A, B)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{
%%     ~~Describe \code{D} here~~
}
  \item{h}{
%%     ~~Describe \code{h} here~~
}
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
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
function (D, h, A, B) 
{
    D[1][which(is.na(D[1]) == TRUE)] <- -Inf
    D[2][which(is.na(D[2]) == TRUE)] <- Inf
    if (D[1] == -Inf && D[2] == Inf) {
        if (is.numeric(A) & is.numeric(B)) {
            a <- A
            b <- B
        }
        else {
            a <- -4
            b <- 4
        }
        ap = diag(attr(numericDeriv(quote(h(a)), "a"), "gradient"))
        bp = diag(attr(numericDeriv(quote(h(b)), "b"), "gradient"))
        if (ap > 0 & bp < 0) {
            T <- c(a, b)
        }
        else if (ap > 0 & bp >= 0) {
            stop("No points to the right of the mode")
        }
        else if (ap <= 0 & bp < 0) {
            stop("No points to the left of the mode")
        }
        else {
            stop("Please give a valid domain for g(x)!\nPlease try again!")
        }
    }
    else if (D[1] == -Inf) {
        a <- -4
        b <- D[2]
        if (a >= b) {
            a <- 2 * b
        }
        T <- c(a, b)
    }
    else if (D[2] == Inf) {
        a <- D[1]
        b <- 4
        if (a >= b) {
            b <- 2 * a
        }
        T <- c(a, b)
    }
    else {
        if (is.numeric(A) & is.numeric(B)) {
            a = A
            b = B
        }
        else {
            a <- D[1] + sqrt(.Machine$double.eps)
            b <- D[2] - sqrt(.Machine$double.eps)
        }
        T <- c(a, b)
    }
    return(T)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
