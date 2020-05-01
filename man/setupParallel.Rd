\name{setupParallel}
\alias{setupParallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setting up parallel computing
}
\description{
This function sets up parallel computing by the package \code{snowfall}.
}
\usage{
setupParallel(ncpus = 1, sourcefile = NULL, sourcelibraries = c("multiridge","penalized","survival","pROC","risksetROC"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ncpus}{
Integer. Number of cpus to use.
}
  \item{sourcefile}{
Character. Additional source files to be loaded in parallel. Only required when parallel computing is also desired for functions
not available in \code{multiridge}.
}
  \item{sourcelibraries}{
Character vector. Libraries to be loaded in parallel. Defaults to the libraries multiridge depends on.
}
}
\details{
Parallel computing is available for several functions that rely on cross-validation. If double CV is used, parallel computing is applied to the outer loop, to optimize efficiency.
}
\value{
No output.
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
Snowfall package for further documentation on parallel computing.
}
\examples{
}
