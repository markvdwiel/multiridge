\name{Scoring}
\alias{Scoring}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate predictions
}
\description{
Evaluates predictions by a score suitable for the corresponding response
}
\usage{
Scoring(lp, Y, model = NULL, score = ifelse(model == "linear", "mse", "loglik"),
  print = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lp}{
Numerical vector. Linear predictor.
}
  \item{Y}{
Response vector: numeric, binary, factor or \code{survival}.
}
  \item{score}{
Character. See Details.
}
  \item{model}{
Character. Any of \code{c("linear", "logistic", "cox")}. Is inferred from
\code{Y} when \code{NULL}.
}
  \item{print}{
Boolean. Should the score be printed on screen.
}
}
\details{
Several scores are allowed, depending on the type of output. For \code{model = "linear"},
\code{score} equals any of \code{c("loglik","mse","abserror","cor","kendall","spearman")}, denoting
CV-ed log-likelihood, mean-squared error, mean absolute error, Pearson (Kendall, Spearman) correlation with response.
For \code{model = "logistic"}, \code{score} equals any of \code{c("loglik","auc", "brier")}, denoting
CV-ed log-likelihood, area-under-the-ROC-curve, and brier score a.k.a. MSE.
For \code{model = "cox"}, \code{score} equals any of \code{c("loglik","cindex")}, denoting
CV-ed log-likelihood, and c-index. %%  ~~ If necessary, more details than the description above ~~
}
\value{
Numerical value.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CVscore}} for obtaining the cross-validated score (for given penalties), and \code{\link{doubleCV}} to obtain doubly cross-validated linear predictors to which \code{Scoring} can be applied to estimated predictive performance by double cross-validation. A full demo and data are available from:\cr
\url{https://drive.google.com/open?id=1NUfeOtN8-KZ8A2HZzveG506nBwgW64e4}
}
