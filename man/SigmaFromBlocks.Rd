\name{SigmaFromBlocks}
\alias{SigmaFromBlocks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create penalized sample cross-product matrix
}
\description{
Creates penalized sample cross-product matrix, dimension \code{nxn}.
}
\usage{
SigmaFromBlocks(XXblocks, penalties, pairing = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XXblocks}{
List of \code{nxn} matrices. Usually output of \code{\link{createXXblocks}}.
}
  \item{penalties}{
Numeric vector, representing penaltyparameters. 
}
  \item{pairing}{
Numerical vector of length 3 or \code{NULL} when pairs are absent. Represents the indices (in \code{XXblocks}) of the two data blocks involved in pairing,
plus the index of the paired block.
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
#Example
#Simulate
Xbl1 <- matrix(rnorm(1000),nrow=10)
Xbl2 <- matrix(rnorm(2000),nrow=10)

#check whether dimensions are correct
ncol(Xbl1)==nrow(Xbl2)

#create cross-product
XXbl <- createXXblocks(list(Xbl1,Xbl2))

#suppose penalties for two data types equal 5,10, respectively
Sigma <- SigmaFromBlocks(XXbl,c(5,10))

#check dimensions (should be n x n)
dim(Sigma)
}
