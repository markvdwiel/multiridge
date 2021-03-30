\name{createXblocks}
\alias{createXblocks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create list of paired data blocks}
\description{
Create list of paired data blocks
}
\usage{
createXblocks(datablocks, which2pair = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datablocks}{
List of data frames or matrices representing \code{b=1,...,B} data blocks of dimensions \code{n x p_b}.
}
  \item{which2pair}{
Integer vector of size 2 (or \code{NULL})
}
}
\details{
Only use this function when you wish to pair two data blocks. If \code{which2pair = NULL} the output
matches the input. If not, the function adds a paired data block, pairing the two data blocks corresponding to the elements of
\code{which2pair}.
}
\value{
List. Same as \code{datablocks} when \code{which2pair = NULL}, or augmented with one paired data block.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{createXXblocks}}
}


