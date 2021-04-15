\name{IWLSridge}
\alias{IWLSridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Iterative weighted least squares algorithm for linear and logistic ridge regression.

}
\description{
Iterative weighted least squares algorithm for linear and logistic ridge regression. Updates the weights and linear predictors until convergence.
}
\usage{
IWLSridge(XXT, Y, X1 = NULL, intercept = TRUE, frac1 = NULL, eps = 1e-07,
maxItr = 25, trace = FALSE, model = NULL, E0 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XXT}{
Matrix. Dimensions \code{nxn}. Sample cross-product from penalized variables, usually computed by \code{\link{SigmaFromBlocks}}.
}
  \item{Y}{
Response vector: numeric, binary, or two-class factor
}
   \item{X1}{
Matrix. Dimension \code{n x p_0, p_0 < n}, representing unpenalized covariates.
}
   \item{intercept}{
Boolean. Should an intercept be included?
}
  \item{frac1}{
Scalar. Prior fraction of cases. Only relevant for \code{model="logistic"}.
}
  \item{eps}{
Scalar. Numerical bound for IWLS convergence.
}
  \item{maxItr}{
Integer. Maximum number of iterations used in IWLS.
}
  \item{trace}{
Boolean. Should the output of the IWLS algorithm be traced?
}
  \item{model}{
Character. Any of \code{c("linear", "logistic")}. Is inferred from
\code{Y} when \code{NULL}. Note that the cox model for survival is covered by the function \code{\link{IWLSCoxridge}}.
}
 \item{E0}{
Numerical vector or \code{NULL}. Optional initial values for linear predictor. Same length as \code{Y}. Usually \code{NULL}, which initializes linear predictor with 0.
}
}
\details{
An (unpenalized) intercept is included by default. To keep the function computationally efficient it returns the linear predictors (which suffice for predictions), instead of parameter estimates. These may be obtained by applying the \code{\link{betasout}} function to the output of this function.}
\value{
List, containing:

\item{etas}{Numerical vector: Final linear predictors}
\item{Ypred}{Predicted survival}
\item{convergence}{Boolean: has IWLS converged?}
\item{nIt}{Number of iterations}
\item{Hres}{Auxiliary list object. Passed on to other functions}
\item{linearized}{Linearized predictions}
\item{unpen}{Boolean: are there any unpenalized covariates involved? Passed on to other functions}
\item{intercept}{Boolean: Is an intercept included? }
\item{eta0}{Numerical vector: Initial linear predictors}
\item{X1}{Matrix: design matrix unpenalized variables}
}

\references{
Mark A. van de Wiel, Mirrelijn van Nee, Armin Rauschenberger (2021). Fast cross-validation for high-dimensional ridge regression. J Comp Graph Stat
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{IWLSCoxridge}} for Cox ridge. \code{\link{betasout}} for obtaining parameter estimates.
\code{\link{predictIWLS}} for predictions on new samples. A full demo and data are available from:\cr
\url{https://drive.google.com/open?id=1NUfeOtN8-KZ8A2HZzveG506nBwgW64e4}
}

\examples{
data(dataXXmirmeth)
resp <- dataXXmirmeth[[1]]
XXmirmeth <- dataXXmirmeth[[2]]
lambdas <- c(100,1000)

# Prepare fitting for the specified penalties.
XXT <- SigmaFromBlocks(XXmirmeth,penalties=lambdas)

# Fit. fit$etas contains the n linear predictors
fit <- IWLSridge(XXT,Y=resp)
}
