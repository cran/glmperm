\name{prr.test}
\alias{prr.test}
\title{
Inference in Generalized Linear Models and the Cox Proportional Hazard Model
}
\description{
A permutation of regressor residuals (prr) test which is used for inference in Generalized Linear Models and Cox Proportional Hazard Models. 
The null hypothesis that the regression coefficient of the covariate of interest is zero is tested.
The procedure is useful when parameter estimates fail to converge or are unreliable due to
small sample size.
}
\usage{
prr.test(formula, var, family=gaussian, data, nrep = 1000, seed=12345, Silent=TRUE, weights,  subset, na.action, start = NULL, etastart, mustart, offset, control = glm.control(...), model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, ...)}
\arguments{
  \item{formula}{a symbolic description of the model to be fitted. If family="Cox" the response variable has to be a survival object 
}
  \item{var}{
the independent variable about which inference
 is to be made; a character string identical to the one in the formular
}
  \item{family}{
a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See family for details of family functions.\\
if family="Cox" the Cox proportional hazard model is fitted (method=coxph.fit).)
}
  \item{data}{
an optional data frame, list or environment containing the variables in the model.
}
  \item{nrep}{
number of Monte Carlo replicates to be used for the permutation test
}
  \item{seed}{
used seed. Default is seed=12345.
}
\item{Silent}{Default is TRUE.}
\item{weights}{see \code{\link{glm}}}
\item{subset}{see \code{\link{glm}}}
\item{na.action}{see \code{\link{glm}}}
\item{start}{see \code{\link{glm}}}
\item{etastart}{see \code{\link{glm}}}
\item{mustart}{see \code{\link{glm}}}
\item{offset}{see \code{\link{glm}}}
\item{control}{see \code{\link{glm}}}
\item{model}{see \code{\link{glm}}}
\item{method}{see \code{\link{glm}}}
\item{x}{see \code{\link{glm}}}
\item{y}{see \code{\link{glm}}}
\item{contrasts}{see \code{\link{glm}}}
\item{...}{see \code{\link{glm}}}
}
\details{
The prr test replaces the independent variable of interest with its residuals from a linear regression on the other independent variables. 
The test statistic for the permutation test is the p-value based on a likelihood ratio test for a generalized linear model or the Cox proportional hazard model. 
Thus, the permutation p-value is the fraction of the permutations that have a likelihood-based p-value less than or equal to that for 
the unpermuted data. Because the p-values for permutations will be a discrete set and round-off errors may occur, prr.test also investigates the sensitivity of the permutation p-value to small variations in the unpermuted likelihood-based p-value.

Missing values are allowed, and observations with them are eliminated.

prr.test calls glm.perm.
}
\value{
\item{fit1}{glm.fit of the full model.}
\item{fit2}{glm.fit of the model without variable 'var'.}
\item{call}{the original function call.}
\item{formula}{the full model formula.}
\item{seed}{used in set.seed(seed).}
\item{Silent}{if Silent=False then the execution time is determined.}
\item{fit1deviance}{the deviance of the full model given in fit1. This value is used for calculating 'LRstat'.}
\item{fit2deviance}{the deviance of the model without the variable 'var' given in 'fit2'. This value is used for calculating 'LRstat'.}
\item{Dispersion}{dispersion factor of the full model. For family binomial and poisson this is set to be 1. 
This value is used for calculating 'LRstat'.
The estimated dispersion factor 'estimated.Dispersion' may differ serverely. 
Then the use of family quasibinomial and quasipoisson might be more appropriate.}
\item{estimated.Dispersion}{the estimated dispersion factor based on the data.}
\item{LRstat}{observed Likelihood Ratio test statistics based on the original data (LRstat=abs(fit1deviance-fit2deviance)/Dispersion).}
\item{p.value.obs}{observed p-value of the Likelihood Ratio test statistcs 'LRstat' calculated via a Chi-square distribution with one degree of freedom.}
\item{p.value.perm}{permutation p-value, fraction of permutation based p-values
 less than or equal to the observed p-value.}
\item{p.value.perm.se}{estimated standard errors associated with permutation p-values}
\item{nobs}{number of observations used}
\item{var}{the independent variable which is tested for association with the response variable.}
\item{nrep.true}{if a fit based on the resampled data failed to converge, a missing value for the corresponding permutation based p-value is produced. 'nrep.true' provides the actual number of permutation based p-values, normally equal or close to 'nrep'.}
}
\references{
Potter D.M. (2005) \emph{A permutation test for inference in logistic regression with small- and moderate-sized datasets}. 
Statistics in Medicine, 24:693-708.
}
\author{
Douglas M. Potter (original version for logistic regression models), Wiebke Werft (extension for generalized linear models and Cox proportional hazard model)
}

\examples{
########### EXAMPLE 1 ##############
## 40 observations, 3 independent variables, family=binomial()
nobs<-40
x1<-rnorm(nobs)
x2<-rnorm(nobs)
xx<-cbind(x1,x2)
x0<-rnorm(nobs)+x1+x2
y<-x0+x1+x2+2*rnorm(nobs)
y<-ifelse(y>0,1,0)
LRtest1 <- prr.test(y~x0+x1+x2,var="x0", family=binomial())
summary(LRtest1)


########### EXAMPLE 2 ##############
## treepipit data from library(coin) for family=poisson()
utils::data(treepipit, package="coin")
LRtest2 <- prr.test(counts ~ cbpiles + coverstorey, var="cbpiles", data=treepipit, family=poisson())
summary(LRtest2)


########### EXAMPLE 3 ##############
## the same as example 2 but with family=quasipoisson()
utils::data(treepipit, package="coin")
LRtest3 <- prr.test(counts ~ cbpiles + coverstorey, var="cbpiles", data=treepipit, family=quasipoisson())
summary(LRtest3)
}
