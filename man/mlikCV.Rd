\name{mlikCV}
\alias{mlikCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Outer-loop cross-validation for estimating performance of marginal likelihood based \code{multiridge}
}
\description{
Outer-loop cross-validation for estimating performance of marginal likelihood based \code{multiridge}.
Outer fold is for testing;  penalty parameter tuning is performed by marginal likelihood estimation}
\usage{
mlikCV(penaltiesinit, XXblocks, Y, pairing = NULL, outfold = 5, nrepeatout = 1,
balance = TRUE,fixedfolds = TRUE,  model = NULL, intercept = ifelse(class(Y) ==
"Surv", FALSE, TRUE),reltol = 1e-04,  trace = FALSE,  optmethod1 = "SANN",
optmethod2 = ifelse(length(penaltiesinit) == 1, "Brent", "Nelder-Mead"),
maxItropt1 = 10, maxItropt2 = 25, parallel = FALSE, pref = NULL,
fixedpen = NULL, sigmasq = 1, opt.sigma=ifelse(model=="linear",TRUE, FALSE))
}


%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{penaltiesinit}{
Numeric vector. Initial values for penaltyparameters. May be obtained from \code{\link{fastCV2}}.
}
  \item{XXblocks}{
List of \code{nxn} matrices. Usually output of \code{\link{createXXblocks}}.
}
  \item{Y}{
Response vector: numeric, binary, factor or \code{survival}.
}
  \item{pairing}{
Numerical vector of length 3 or \code{NULL} when pairs are absent. Represents the indices (in \code{XXblocks}) of the two data blocks involved in pairing,
plus the index of the paired block.
}
  \item{outfold}{
Integer. Outer fold for test samples.
}
  \item{nrepeatout}{
Integer. Number of repeated splits for outer fold.
}
  \item{balance}{
Boolean. Should the splits be balanced in terms of response labels?
}
  \item{fixedfolds}{
Boolean. Should fixed splits be used for reproducibility?
}
   \item{intercept}{
Boolean. Should an intercept be included?
}
  \item{model}{
Character. Any of \code{c("linear", "logistic", "cox")}. Is inferred from
\code{Y} when \code{NULL}.
}
  \item{trace}{
Boolean. Should the output of the IWLS algorithm be traced?
}
  \item{reltol}{
Scalar. Relative tolerance for optimization methods.
}
  \item{optmethod1}{
Character. First, global search method. Any of the methods \code{c("Brent", "Nelder-Mead", "Sann")} may be used, but
simulated annealing by \code{"Sann"} is recommended to search a wide landscape. Other unconstrained methods
offered by \code{\link{optim}} may also be used, but have not been tested.
}
  \item{optmethod2}{
Character. Second, local search method. Any of the methods \code{c("Brent", "Nelder-Mead", "Sann")} may be used, but
\code{"Nelder-Mead"} is generally recommended. Other unconstrained methods
offered by \code{\link{optim}} may also be used, but have not been tested.
}
  \item{maxItropt1}{
Integer. Maximum number of iterations for \code{optmethod1}.
}
  \item{maxItropt2}{
Integer. Maximum number of iterations for \code{optmethod2}.
}
  \item{parallel}{
Boolean. Should computation be done in parallel? If \code{TRUE}, requires to run \code{\link{setupParallel}} first.
}
  \item{pref}{
Integer vector or \code{NULL}. Contains indices of data types in \code{XXblocks} that are preferential.
}
  \item{fixedpen}{
Integer vector or \code{NULL}. Contains indices of data types of which penalty is fixed to the corresponding value in \code{penaltiesinit}.
}
\item{sigmasq}{Default error variance.
}
\item{opt.sigma}{Boolean. Should the error variance be optimized as well? Only relevant for \code{model="linear"}.
}


}
\details{
WARNING: this function may be very time-consuming. The number of evaluations may equal \code{nrepeatout*outerfold*(maxItropt1+maxItropt2)}. Computing time may be estimated by multiplying computing time of \code{\link{optLambdas_mgcvWrap}} by
\code{nrepeatout*outerfold}.
}
\value{
List with the following components:
\item{sampleindex}{Numerical vector: sample indices}
\item{true}{True responses}
\item{linpred}{Cross-validated linear predictors}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{optLambdas_mgcv}}, \code{\link{optLambdas_mgcvWrap}} which optimize the penalties.
\code{\link{Scoring}} which may applied to output of this function to obtain overall cross-validated performance score.
\code{\link{doubleCV}} for double cross-validation counterpart. A full demo and data are available from:\cr
\url{https://drive.google.com/open?id=1NUfeOtN8-KZ8A2HZzveG506nBwgW64e4}
}

\examples{
data(dataXXmirmeth)
resp <- dataXXmirmeth[[1]]
XXmirmeth <- dataXXmirmeth[[2]]

# Find initial lambdas: fast CV per data block separately.
cvperblock2 <- fastCV2(XXblocks=XXmirmeth,Y=resp,kfold=10,fixedfolds = TRUE)
lambdas <- cvperblock2$lambdas

# Outer cross-validation, inner marginal likelihood optimization
\dontrun{
perfmlik <- mlikCV(penaltiesinit=lambdas,XXblocks=XXmirmeth,Y=resp,outfold=10,
nrepeatout=1)


# Performance metrics
Scoring(perfmlik$linpred,perfmlik$true,score="auc",print=TRUE)
Scoring(perfmlik$linpred,perfmlik$true,score="brier",print=TRUE)
Scoring(perfmlik$linpred,perfmlik$true,score="loglik",print=TRUE)
}
}

