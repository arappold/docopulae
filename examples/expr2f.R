library(copula)

## build f
margins = list(alist(pdf=dnorm(y1, mu1, sd1),
                     cdf=pnorm(y1, mu1, sd1)),
               alist(pdf=dnorm(y2, mu2, sd2),
                     cdf=pnorm(y2, mu2, sd2)))
C = claytonCopula()
f = buildf(margins, C)
f

## turn f into function
expr2f(f, yMap=list(y1=1, y2=2),
       thetaMap=list(alpha='alpha', mu1='mu1', mu2='mu2',
                                    sd1='sd1', sd2='sd2'))
