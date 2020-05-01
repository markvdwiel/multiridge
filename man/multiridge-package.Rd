\name{multiridge-package}
\alias{multiridge-package}
\alias{multiridge}
\docType{package}
\title{
Fast cross-validation for multi-penalty ridge regression
}
\description{
The package implements multi-penalty linear, logistic and cox ridge regression, including estimation of the penalty parameters by efficient (repeated) cross-validation. It allows for multiple high-dimensional data types that require penalization, as well as unpenalized variables. Moreover, it allows a paired penalty for paired data types, and preferential data types can be specified. 
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
\code{\link{fastCV}}: Fast cross-validation per data block\cr
\code{\link{IWLSCoxridge}}: Iterative weighted least squares algorithm for Cox ridge regression \cr
\code{\link{IWLSridge}}: Iterative weighted least squares algorithm for linear and logistic ridge regression \cr
\code{\link{optLambdasWrap}}: Find optimal ridge penalties  \cr
\code{\link{predictIWLS}}: Predictions from ridge fits \cr
\code{\link{setupParallel}}: Setting up parallel computing\cr
\code{\link{SigmaFromBlocks}}: Create penalized sample cross-product matrix
}

\author{
Mark A. van de Wiel (mark.vdwiel@amsterdamumc.nl)
}
\references{
Mark A. van de Wiel, Mirrelijn van Nee, Armin Rauschenberger (2020). Fast CV for multi-omics ridge regression.
}
\keyword{ package }
\seealso{
}
\examples{
#A complete demo on real data is available at: 
}
