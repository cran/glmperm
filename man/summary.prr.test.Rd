\name{summary.prr.test}
\alias{summary.prr.test}
\title{
Summary of the prr.test results.
}
\description{
summary is a generic function used to produce result summaries of objects returned by prr.test
}
\usage{
\method{summary}{prr.test}(object, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{object}{
a result object of the call prr.test.
}
  \item{digits}{
how many digits the values of the results should have.
}
\item{...}{
... additional arguments affecting the summary produced. 
}
}
\details{
summary.prr.test prints an overview of the main results of the permutation Likelihood Ratio Test of the prr.test. 
}

\author{
Wiebke Werft 
}

\examples{
utils::data(anorexia, package="MASS")
LRtest <- prr.test(Postwt ~ Prewt + Treat, var="Prewt", data = anorexia, family = gaussian())
summary(LRtest)
  }


