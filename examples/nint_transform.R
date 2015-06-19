library(docopulae)

f = sin
s = nint_space(nint_intvDim(-2.11, 1.76))
nint_integrate(f, s)

tt = nint_transform(f, s, 1, 'tan')
tt$space
nint_integrate(tt$f, tt$space)

tt = nint_transform(f, s, 1, 'ratio')
tt$space
nint_integrate(tt$f, tt$space)

# probability integral transform
tt = nint_transform(f, s, 1, list(g=pnorm,
    gij=function(x) { t1 = qnorm(x); cbind(t1, 1/dnorm(t1)) }))
tt$space
nint_integrate(tt$f, tt$space)


f = function(x) prod(1/(1 + x**2))
s = nint_space(nint_intvDim(-1, Inf),
               nint_intvDim(-Inf, Inf))
s

nint_integrate(f, s) # stats::integrate takes care of Inf limits

tt = nint_transform(f, s, 1:2, 'tan')
tt$space
nint_integrate(tt$f, tt$space)

tt = nint_transform(f, s, 1:2, 'ratio')
tt$space
nint_integrate(tt$f, tt$space)

# probability integral transform
tt = nint_transform(f, s, 1:2, list(g=pnorm,
    gij=function(x) { t1 = qnorm(x); cbind(t1, 1/dnorm(t1)) }))
tt$space
# won't terminate. I guess the transformation is too steep
#nint_integrate(tt$f, tt$space)

