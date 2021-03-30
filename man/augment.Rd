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




