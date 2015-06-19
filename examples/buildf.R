library(docopulae)
library(copula)
library(mvtnorm)

## build bivariate normal pdf
margins = function(y, theta) {
    mu = c(theta$mu1, theta$mu2)
    cbind(dnorm(y, mu), pnorm(y, mu))
}
C = copula::normalCopula()

f = buildf(margins, C, 'alpha')

## plot density
theta = list(mu1=2, mu2=-3, alpha=0.4)
y1 = seq(0, 4, length.out=51)
y2 = seq(-5, -1, length.out=51)
z = outer(y1, y2, function(y1, y2) apply(cbind(y1, y2), 1, f, theta))
contour(y1, y2, z)

## add theoretical density
C@parameters = theta$alpha
zz = outer(y1, y2,
    function(y1, y2) dmvnorm(cbind(y1 - theta$mu1, y2 - theta$mu2),
                             sigma=getSigma(C)))
contour(y1, y2, zz, col='red', add=TRUE)

## build bivariate pdf with normal margins and clayton copula
margins = list(alist(pdf=dnorm(y1, mu1, 1),
                     cdf=pnorm(y1, mu1, 1)),
               alist(pdf=dnorm(y2, mu2, 1),
                     cdf=pnorm(y2, mu2, 1)))
C = claytonCopula()
f1 = buildf(margins, C)
f1
ff1 = expr2f(f1, yMap=list(y1=1, y2=2),
                 thetaMap=list(mu1='mu1', mu2='mu2', alpha='alpha'))
ff1

f2 = buildf(function(y, theta) {
    mu = c(theta$mu1, theta$mu2)
    cbind(dnorm(y, mu, 1), pnorm(y, mu, 1))
}, C, 'alpha')

## plot both densities
theta = list(mu1=2, mu2=-3, alpha=2) # tau = 0.5

y1 = seq(0, 4, length.out=51)
y2 = seq(-5, -1, length.out=51)
z = outer(y1, y2, function(y1, y2) apply(cbind(y1, y2), 1, ff1, theta))
contour(y1, y2, z)

zz = outer(y1, y2, function(y1, y2) apply(cbind(y1, y2), 1, f2, theta))
contour(y1, y2, zz, col='red', add=TRUE)
