\name{CVfolds}
\alias{CVfolds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates (repeated) cross-validation folds
}
\description{
Creates (repeated) cross-validation folds for samples
}
\usage{
CVfolds(Y, model = NULL, balance = TRUE, kfold = 10, fixedfolds = TRUE, nrepeat = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
Response vector: numeric, binary, factor or \code{survival}.
}
  \item{model}{
Character. Any of \code{c("linear", "logistic", "cox")}. Is inferred from
\code{Y} when \code{NULL}.
}
  \item{balance}{
Boolean. Should the splits be balanced in terms of response labels?
}
  \item{kfold}{
Integer. Desired fold.
}
  \item{fixedfolds}{
Boolean. Should fixed splits be used for reproducibility?
}
  \item{nrepeat}{
Numeric. Number of repeats.
}
}
\details{
Creates (repeated), possibly balanced,  splits of the samples. Computing time will often largely depend on
on \code{kfold*nrepeat}, the number of training-test splits evaluated.
}
\value{
List object with \code{kfold*nrepeat} elements containing the sample indices of the left-out samples per split.
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
