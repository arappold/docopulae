library(docopulae)
library(copula)

eta1 = quote(beta1 + beta2*x + beta3*x**2)
eta2 = quote(beta4*x + beta5*x**3 + beta6*x**4)
margins = list(list(pdf=substitute(dnorm(y1, eta1, 1), list(eta1=eta1)),
                    cdf=substitute(pnorm(y1, eta1, 1), list(eta1=eta1))),
               list(pdf=substitute(dnorm(y2, eta2, 1), list(eta2=eta2)),
                    cdf=substitute(pnorm(y2, eta2, 1), list(eta2=eta2))))
C = frankCopula()
f = buildf(margins, C)
f

n = c('alpha', 'beta1', 'beta2', 'beta3', 'beta4', 'beta5', 'beta6')
d2logf = Deriv2Logf(f, n)
str(d2logf)



