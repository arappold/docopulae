library(docopulae)
library(copula)

## build f
margins = list(alist(pdf=dnorm(y1, beta10 + beta11*x, 1),
                     cdf=pnorm(y1, beta10 + beta11*x, 1)),
               alist(pdf=dnorm(y2, beta20 + beta21*x, 1),
                     cdf=pnorm(y2, beta20 + beta21*x, 1)))
C = claytonCopula()
f = buildf(margins, C)
f

theta = list(alpha=2, beta10=1, beta11=2, beta20=3, beta21=4, x=5)

## build first derivative
n = c('alpha', 'beta10', 'beta11', 'beta20', 'beta21')
yMap = list(y1=1, y2=2)
thetaMap = list(alpha='alpha', beta10='beta10', beta11='beta11',
                               beta20='beta20', beta21='beta21', x='x')
dlogf = DerivLogf(f, n, yMap=yMap, thetaMap=thetaMap)
str(dlogf)

## build reference
margins = function(y, theta) {
    mu = c(theta$beta10 + theta$beta11*theta$x,
           theta$beta20 + theta$beta21*theta$x)
    cbind(dnorm(y, mu, 1), pnorm(y, mu, 1))
}
f2 = buildf(margins, C, 'alpha')
f2

dlogf2 = numDerivLogf(f2)

## evaluate
y = c(11, 22)
m1 = quote( sapply(n, function(i) dlogf(y, theta, i)) )
m2 = quote( sapply(n, function(i) dlogf2(y, theta, i)) )
eval(m1)
eval(m2)

system.time(replicate(50, eval(m1)))
system.time(replicate(50, eval(m2)))

## build second derivative
d2logf = Deriv2Logf(f, n, yMap=yMap, thetaMap=thetaMap)
str(d2logf)

## build reference
d2logf2 = numDeriv2Logf(f2)

## evaluate
m1 = quote( outer(n,n, function(a, b) apply(cbind(a, b), 1,
    function(x) d2logf(y, theta, x[1], x[2]))) )
m2 = quote( outer(n, n, function(a, b) apply(cbind(a, b), 1,
    function(x) d2logf2(y, theta, x[1], x[2]))) )
eval(m1)
eval(m2)

## benchmark
system.time(replicate(50, eval(m1)))
system.time(replicate(5, eval(m2))) * 10
