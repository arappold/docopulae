library(docopulae)

## integrate with defaults
nint_integrate(sin, nint_space(nint_intvDim(0, 2*pi)))

## replace nint_integrateNCube
f = function(f, lowerLimit, upperLimit, ...) {
    print('using integrateA')
    integrateA(f, lowerLimit, upperLimit, ..., subdivisions=2)
}
f = nint_integrateNCube_integrate(f)
unlockBinding('nint_integrateNCube', environment(nint_integrate))
assign('nint_integrateNCube', f, envir=environment(nint_integrate))

## integrate with replacement
nint_integrate(sin, nint_space(nint_intvDim(0, 2*pi))) # no difference here

## replace nint_integrateNCube
f = function(f, lowerLimit, upperLimit, ...) {
    print('using cubature')
    r = cubature::adaptIntegrate(f, lowerLimit, upperLimit, ..., maxEval=1e3)
    return(r$integral)
}
unlockBinding('nint_integrateNCube', environment(nint_integrate))
assign('nint_integrateNCube', f, envir=environment(nint_integrate))

## integrate with replacement
nint_integrate(sin, nint_space(nint_intvDim(0, 2*pi)))

## replace nint_integrateNCube
f = function(dimension) {
    print('using SparseGrid')
    SparseGrid::createIntegrationGrid('GQU', dimension, 7)
}
f = nint_integrateNCube_SparseGrid(f)
unlockBinding('nint_integrateNCube', environment(nint_integrate))
assign('nint_integrateNCube', f, envir=environment(nint_integrate))

## integrate with replacement
nint_integrate(sin, nint_space(nint_intvDim(0, 2*pi)))
