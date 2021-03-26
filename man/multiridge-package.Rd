\name{multiridge-package}
\alias{multiridge-package}
\alias{multiridge}
\docType{package}
\title{
Fast cross-validation for multi-penalty ridge regression
}
\description{
The package implements multi-penalty linear, logistic and cox ridge regression, including estimation of the penalty parameters by efficient (repeated) cross-validation or marginal likelihood maximization. It allows for multiple high-dimensional data types that require penalization, as well as unpenalized variables. Moreover, it allows a paired penalty for paired data types, and preferential data types can be specified.
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{multiridge}
\packageIndices{multiridge}
\cr
\code{\link{betasout}}: Coefficient estimates from (converged) IWLS fit \cr
\code{\link{createXXblocks}}: Creates list of (unscaled) sample covariance matrices \cr
\code{\link{CVscore}}: Cross-validated score for given penalty parameters \cr
\code{\link{doubleCV}}:  Double cross-validation for estimating performance \cr
\code{\link{fastCV}}: Fast, SVD-based cross-validation per data block; based on \code{penalized} package\cr
\code{\link{fastCV2}}: Fast cross-validation per data block; no dependency\cr
\code{\link{IWLSCoxridge}}: Iterative weighted least squares algorithm for Cox ridge regression \cr
\code{\link{IWLSridge}}: Iterative weighted least squares algorithm for linear and logistic ridge regression \cr
\code{\link{mlikCV}}:  Cross-validation for estimating performance of marginal likelihood estimation\cr
\code{\link{optLambdasWrap}}: Find optimal ridge penalties by cross-validation  \cr
\code{\link{optLambdas_mgcvWrap}}: Find optimal ridge penalties in terms of marginal likelihood  \cr
\code{\link{predictIWLS}}: Predictions from ridge fits \cr
\code{\link{setupParallel}}: Setting up parallel computing\cr
\code{\link{SigmaFromBlocks}}: Create penalized sample cross-product matrix
}

\author{
Mark A. van de Wiel (mark.vdwiel@amsterdamumc.nl)
}
\references{
Mark A. van de Wiel, Mirrelijn van Nee, Armin Rauschenberger (2021). Fast cross-validation for high-dimensional ridge regression. J Comp Graph Stat
}
\keyword{ package }
\seealso{A full demo and data are available from:
\url{https://drive.google.com/open?id=1NUfeOtN8-KZ8A2HZzveG506nBwgW64e4}
}
\examples{
}
