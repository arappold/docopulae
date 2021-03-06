% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{buildf}
\alias{buildf}
\title{Build probability density or mass Function}
\usage{
buildf(margins, continuous, copula, parNames = NULL,
  simplifyAndCache = T)
}
\arguments{
\item{margins}{either \itemize{
\item \code{function(y, theta, ...)}, where \code{theta} is a list of parameters.
It shall return a column matrix of two, the probability densities and cumulative distributions.
\item a list of pairs of expressions, each named \code{"pdf"} and \code{"cdf"}, the probability density and cumulative distribution.
}}

\item{continuous}{\code{TRUE} if margins are continuous. See details.}

\item{copula}{if \code{margins} is \itemize{
\item a function then either a copula object from package \pkg{copula} or \code{function(u, theta, ...)}, a probability density function if \code{continuous} else a cumulative distribution function.
\item a list then either a copula object from package \pkg{copula} which contains distribution expressions or an expression for the probability density if \code{continuous} else the cumulative distribution which uses \code{u1},\code{u2},...
}}

\item{parNames}{if \itemize{
\item (optional) \code{margins} is a function and \code{copula} is a copula object then a vector of names or indices, the sequence of copula parameters in \code{theta}.
\code{0} or \code{""} identifies copula parameters to skip.
\item \code{margins} is a list and \code{copula} is a copula object then a named list of names or indices, mapping parameters in \code{theta} to copula parameter variables.
See \code{copula@exprdist}.
}}

\item{simplifyAndCache}{(if \code{margins} is a list) simplify and cache the result using \code{\link[Deriv]{Simplify}} and \code{\link[Deriv]{Cache}} from package \pkg{Deriv} if available.}
}
\value{
\code{buildf} returns \code{function(y, theta, ...)}, the joint probability density or mass function.
}
\description{
\code{buildf} builds a joint probability density or mass function from marginal distributions and a copula.
}
\details{
Please note that expressions are not validated.

If \code{continuous} is \code{FALSE}, dimensionality shall be 2 and both dimensions shall be discrete.
The joint probability mass is defined by
\deqn{C(F_{1}(y_{1}),F_{2}(y_{2}))-C(F_{1}(y_{1}-1),F_{2}(y_{2}))-C(F_{1}(y_{1}),F_{2}(y_{2}-1))+C(F_{1}(y_{1}-1),F_{2}(y_{2}-1))}{C(F1(y[1]), F2(y[2])) - C(F1(y[1] - 1), F2(y[2])) - C(F1(y[1]), F2(y[2] - 1)) + C(F1(y[1] - 1), F2(y[2] - 1))}
where \eqn{C}, \eqn{F_{1}}{F1}, and \eqn{F_{2}}{F2} depend on \eqn{\theta} and \eqn{y_{i}\ge0}{y[i] \ge 0}.
}
\examples{
\donttest{
## for an actual use case see examples for param

library(copula)
library(mvtnorm)

## build bivariate normal
margins = function(y, theta) {
    mu = c(theta$mu1, theta$mu2)
    cbind(dnorm(y, mean=mu, sd=1), pnorm(y, mean=mu, sd=1))
}
copula = normalCopula()

# args: function, copula object, parNames
f1 = buildf(margins, TRUE, copula, parNames='alpha1')
f1 # uses theta[['alpha1']] as copula parameter

## evaluate and plot
theta = list(mu1=2, mu2=-3, alpha1=0.4)

y1 = seq(0, 4, length.out=51)
y2 = seq(-5, -1, length.out=51)
v1 = outer(y1, y2, function(z1, z2) apply(cbind(z1, z2), 1, f1, theta))
str(v1)
contour(y1, y2, v1, main='f1', xlab='y1', ylab='y2')

## compare with bivariate normal from mvtnorm
copula@parameters = theta$alpha1
v = outer(y1, y2, function(yy1, yy2)
    dmvnorm(cbind(yy1, yy2), mean=c(theta$mu1, theta$mu2),
                             sigma=getSigma(copula)))
all.equal(v1, v)


## build bivariate pdf with normal margins and Clayton copula
margins = list(list(pdf=quote(dnorm(y[1], theta$mu1, 1)),
                    cdf=quote(pnorm(y[1], theta$mu1, 1))),
               list(pdf=quote(dnorm(y[2], theta$mu2, 1)),
                    cdf=quote(pnorm(y[2], theta$mu2, 1))))
copula = claytonCopula()

# args: list, copula object, parNames
f2 = buildf(margins, TRUE, copula, list(alpha='alpha1'))
f2

## evaluate and plot
theta = list(mu1=2, mu2=-3, alpha1=2)

y1 = seq(0, 4, length.out=51)
y2 = seq(-5, -1, length.out=51)
v2 = outer(y1, y2, function(z1, z2) apply(cbind(z1, z2), 1, f2, theta))
str(v2)
contour(y1, y2, v2, main='f2', xlab='y1', ylab='y2')

## build alternatives
cexpr = substituteDirect(copula@exprdist$pdf,
                         list(alpha=quote(theta$alpha1)))
# args: list, expression
f3 = buildf(margins, TRUE, cexpr) # equivalent to f2
f3

margins = function(y, theta) {
    mu = c(theta$mu1, theta$mu2)
    cbind(dnorm(y, mean=mu, sd=1), pnorm(y, mean=mu, sd=1))
}
# args: function, copula object, parNames
f4 = buildf(margins, TRUE, copula, 'alpha1')
f4

cpdf = function(u, theta) {
    copula@parameters = theta$alpha1
    dCopula(u, copula)
}
# args: function, function
f5 = buildf(margins, TRUE, cpdf) # equivalent to f4
f5

# args: function, copula object
copula@parameters = 2
f6 = buildf(margins, TRUE, copula)
f6 # uses copula@parameters

cpdf = function(u, theta) dCopula(u, copula)
# args: function, function
f7 = buildf(margins, TRUE, cpdf) # equivalent to f6
f7

## compare all
vv = lapply(list(f3, f4, f5, f6, f7), function(f)
    outer(y1, y2, function(z1, z2) apply(cbind(z1, z2), 1, f, theta)))
sapply(vv, all.equal, v2)
}
}
\seealso{
\pkg{copula}, \code{\link[Deriv]{Simplify}}, \code{\link[Deriv]{Cache}}, \code{\link{numDerivLogf}}, \code{\link{DerivLogf}}, \code{\link{fisherI}}
}
