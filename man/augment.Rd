\name{augment}
\alias{augment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Augment data with zeros.
}
\description{
This function augments data with zeros to allow pairing of data on the same variables, but from DIFFERENT samples}
\usage{
augment(Xdata1, Xdata2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xdata1}{
Data frame or data matrix of dimension \code{n_1 x p}.
}
  \item{Xdata2}{
Data frame or data matrix of dimension \code{n_2 x p}
}
  }
\details{
Xdata1 and Xdata2 should have the same number of columns. These columns represent variables. Augments both data matrices with zeros,
such that the matrices can be paired using \code{\link{createXXblocks}} on the output of this function.
}
\value{
List
\item{Xaug1}{Augmented data matrix 1}
\item{Xaug2}{Augmented data matrix 2}
}

\examples{
#Example
#Simulate
n1 <- 10
n2 <- 20
p <- 100
X1 <- matrix(rnorm(p*n1),nrow=n1)
X2 <- matrix(rnorm(p*n2),nrow=n2)

#check whether column dimension is correct
ncol(X1)==ncol(X2)

#create cross-product
Xaugm <- augment(X1,X2)

#check dimensions (should be (n1+n2) x p)
dim(Xaugm[[1]])
dim(Xaugm[[2]])
}




