\name{createXXblocks}
\alias{createXXblocks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates list of (unscaled) sample covariance matrices
}
\description{
Creates list of (unscaled) sample covariance matrices \code{X_b \%*\% t(X_b)} for data blocks b = 1,..., B.}
\usage{
createXXblocks(datablocks, datablocksnew = NULL, which2pair = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datablocks}{
List of data frames or matrices
}
  \item{datablocksnew}{
List of data frames or matrices
}
  \item{which2pair}{
Integer vector of size 2 (or \code{NULL})
}
}
\details{
The efficiency of \code{multiridge} for high-dimendional data relies largely on this function:
all iterative calculation are performed on the out put of this function, which contains \code{B} blocks of
\code{nxn} matrices. If \code{which2pair != NULL}, the function adds a paired covariance block, pairing the two data blocks corresponding to the elements of \code{which2pair}. If predictions for new samples are desired, one also needs to specify
\code{datablocksnew}, which should have he exact same format as \code{datablocks} with matching column dimension (number of variables).
}
\value{
List. Same number of component as \code{datablocks} when \code{which2pair = NULL}, or augmented with one paired data block.
Dimension is \code{nxn} for all components.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{createXblocks}}, which is required when parameter estimates are desired (not needed for prediction). A full demo and data are available from:\cr
\url{https://drive.google.com/open?id=1NUfeOtN8-KZ8A2HZzveG506nBwgW64e4}
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

