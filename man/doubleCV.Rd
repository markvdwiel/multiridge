\name{doubleCV}
\alias{doubleCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Double cross-validation for estimating performance of \code{multiridge}
}
\description{
Double cross-validation for estimating performance of \code{multiridge}. Outer fold is for testing, inner fold for penalty parameter tuning}
\usage{
doubleCV(penaltiesinit, XXblocks, Y, X1 = NULL, pairing = NULL, outfold = 5,
  infold = 10, nrepeatout =   1, nrepeatin = 1, balance = TRUE, fixedfolds =
  TRUE, intercept = ifelse(class(Y) == "Surv", FALSE,     TRUE), frac1 = NULL,
  score = "loglik",model = NULL, eps = 1e-07, maxItr = 10, trace = FALSE,
  printCV   = TRUE, reltol = 1e-04, optmethod1 = "SANN", optmethod2 =
  ifelse(length(penaltiesinit) == 1, "Brent", "Nelder-Mead"), maxItropt1 = 10,
  maxItropt2 = 25, save = FALSE, parallel = FALSE, pref = NULL, fixedpen = NULL)
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
  \item{X1}{
Matrix. Dimension \code{n x p_0, p_0 < n}, representing unpenalized covariates }
  \item{pairing}{
Numerical vector of length 3 or \code{NULL} when pairs are absent. Represents the indices (in \code{XXblocks}) of the two data blocks involved in pairing,
plus the index of the paired block.
}
  \item{outfold}{
Integer. Outer fold for test samples.
}
  \item{infold}{
Integer. Inner fold for tuning penalty parameters.
}
  \item{nrepeatout}{
Integer. Number of repeated splits for outer fold.
}
  \item{nrepeatin}{
Integer. Number of repeated splits for inner fold.
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
  \item{frac1}{
Scalar. Prior fraction of cases. Only relevant for \code{model=" logistic"}.
}
  \item{score}{
Character. See Details.
}
  \item{model}{
Character. Any of \code{c("linear", "logistic", "cox")}. Is inferred from
\code{Y} when \code{NULL}.
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
  \item{printCV}{
Boolean. Should the CV-score be printed on screen?
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
\item{save}{
Boolean. If TRUE appends the penalties and resulting CVscore to global variable \code{allscores}
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
}
\details{
WARNING: this function may be very time-consuming. The number of evaluations may equal \code{nrepeatout*outerfold*nrepeatin*innerfold*maxItr*(maxItropt1+maxItropt2)}. Computing time may be estimated by multiplying computing time of \code{\link{optLambdasWrap}} by
\code{nrepeatout*outerfold}. See \code{\link{Scoring}} for details on \code{score}.
}
\value{
List with the following components:
\item{sampleindex}{Numerical vector: sample indices}
\item{true}{True responses}
\item{linpred}{Cross-validated linear predictors}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{optLambdas}}, \code{\link{optLambdasWrap}} which optimize the penalties.
\code{\link{Scoring}} which may applied to output of this function to obtain overall cross-validated performance score. A full demo and data are available from:\cr
\url{https://drive.google.com/open?id=1NUfeOtN8-KZ8A2HZzveG506nBwgW64e4}
}

\examples{
data(dataXXmirmeth)
resp <- dataXXmirmeth[[1]]
XXmirmeth <- dataXXmirmeth[[2]]

# Find initial lambdas: fast CV per data block separately.
cvperblock2 <- fastCV2(XXblocks=XXmirmeth,Y=resp,kfold=10,fixedfolds = TRUE)
lambdas <- cvperblock2$lambdas

# Double cross-validation
\dontrun{
perf <- doubleCV(penaltiesinit=lambdas,XXblocks=XXmirmeth,Y=resp,
score="loglik",outfold=10, infold=10, nrepeatout=1, nrepeatin=3, parallel=TRUE)

# Performance metrics
Scoring(perf$linpred,perf$true,score="auc",print=TRUE)
Scoring(perf$linpred,perf$true,score="brier",print=TRUE)
Scoring(perf$linpred,perf$true,score="loglik",print=TRUE)
}
}
