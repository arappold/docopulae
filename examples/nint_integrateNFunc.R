library(docopulae)

## volume of sphere
s = nint_space(
    nint_intvDim(-1, 1),
    nint_funcDim(function(x) nint_intvDim(c(-1, 1) * sin(acos(x[1])) )),
    nint_funcDim(function(x) {
        r = sin(acos(x[1]))
        nint_intvDim(c(-1, 1) * r*cos(asin(x[2] / r)))
    }))
nint_integrate(function(x) 1, s) # 4*pi/3

## replace nint_integrateNFunc
f = function(f, funcs, x0, i0, ...) {
    print(list(f=f, funcs=funcs, x0=x0, i0=i0, ...))
    # disregard unknown structure of f
    stop('do something')
}
unlockBinding('nint_integrateNFunc', environment(nint_integrate))
assign('nint_integrateNFunc', f, envir=environment(nint_integrate))

## integrate with replacement
try(nint_integrate(function(x) 1, s))
