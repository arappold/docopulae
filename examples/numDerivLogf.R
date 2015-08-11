library(copula)

## build f
margins = function(y, theta) {
    mu = c(theta$mu1, theta$mu2)
    cbind(dnorm(y, mu, 1), pnorm(y, mu, 1))
}
C = copula::frankCopula()

f = buildf(margins, C, 'alpha')
f

theta = list(alpha=-2, mu1=3, mu2=-8)
y = c(3, -8)
f(y, theta)

## build first derivative
dlogf = numDerivLogf(f)
dlogf

dlogf(y, theta, 'alpha')
n = names(theta)
sapply(n, function(i) dlogf(y, theta, i))

## build second derivative
d2logf = numDeriv2Logf(f)
d2logf

d2logf(y, theta, 'alpha', 'mu1')
m = outer(n, n, function(a, b) apply(cbind(a, b), 1,
    function(x) d2logf(y, theta, x[1], x[2])))
dimnames(m) = list(n, n)
m
