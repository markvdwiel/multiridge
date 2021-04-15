\name{dataXXmirmeth}
\alias{dataXXmirmeth}
\docType{data}
\title{
Contains R-object \code{dataXXmirmeth}
}
\description{This list object contains the binary response (control/case) and two data blocks corresponding to miRNA and methylation data
}
\usage{data(dataXXmirmeth)}
\format{
The format is a list with two components:
resp: numeric (0/1) [1:43]\\
XXmirmeth: list with 2 components, each a matrix [1:43,1:43]\\
}
\details{
The object \code{XXmirmeth} is created by applying \code{createXXblocks(list(datamir,datameth))}, where
objects \code{datamir} and \code{datameth} are large data matrices stored in the
\code{mirmethdata.Rdata} file, which is available from the link below.
}

\source{
Snoek, B. C. et al. (2019), Genome-wide microRNA analysis of HPV-positive self-samples yields novel
triage markers for early detection of cervical cancer, International Journal of Cancer
144(2), 372-379.

Verlaat, W. et al. (2018), Identification and validation of a 3-gene methylation classifier for hpv-based cervical screening on self-samples, Clinical Cancer Research 24(14), 3456-3464.
}

\references{
Mark A. van de Wiel, Mirrelijn van Nee, Armin Rauschenberger (2021). Fast cross-validation for multi-penalty
high-dimensional ridge regression. J Comp Graph Stat
}
\examples{
data(dataXXmirmeth)
resp <- dataXXmirmeth[[1]]
XXmirmeth <- dataXXmirmeth[[2]]
}

\seealso{
\code{createXXblocks}. Source data file is available from:\cr
\url{https://drive.google.com/open?id=1NUfeOtN8-KZ8A2HZzveG506nBwgW64e4}
}

\keyword{datasets}
