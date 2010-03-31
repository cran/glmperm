\name{glmperm-package}
\Rdversion{1.1}
\alias{glmperm-package}
\alias{glmperm}
\docType{package}
\title{
Inference in Generalized Linear Models
}
\description{
A permutation Likelihood Ratio test is used for inference in generalized linear models.
The procedure is useful when parameter estimates fail to converge or are unreliable due to
small sample size. 
}
\details{
\tabular{ll}{
Package: \tab glmperm\cr
Type: \tab Package\cr
Version: \tab 1.0-2\cr
Date: \tab 2009-31-03\cr
License: \tab GPL-3\cr
}
The function prr.test requires that the user identifies one independent 
variable for which a test of the null hypothesis (that the regression
coefficient is zero) is to be performed. This variable of interest is supplied as a character string 'var'; the full model is supplied as a formula.}
\author{
Douglas M. Potter (original version for logistic regression models), Wiebke Werft (extension for generalized linear models)

Maintainer: Wiebke Werft <w.werft@dkfz.de>
}
\references{
Potter D.M. (2005) \emph{A permutation test for inference in 
logistic regression with small- and moderate-sized datasets}. 
Statistics in Medicine, 24:693-708.
}
\examples{
### Example 1 ###
nobs<-40
x1<-rnorm(nobs)
x2<-rnorm(nobs)
xx<-cbind(x1,x2)
x0<-rnorm(nobs)+x1+x2
y<-x0+x1+x2+2*rnorm(nobs)
y<-ifelse(y>0,1,0)
LRtest1 <- prr.test(y~x0+x1+x2,var="x0", family=binomial())
summary(LRtest1)

### Example 2 ###
utils::data(anorexia, package="MASS")
LRtest <- prr.test(Postwt ~ Prewt + Treat, var="Prewt", data = anorexia, family = gaussian())
summary(LRtest)
}
