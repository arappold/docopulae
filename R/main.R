

#' Parametric Model
#'
#' \code{param} creates an initial parametric model object.
#' Unlike other model statements this function does not perform any computation.
#'
#' @param fisherIf \code{function(x, ...)}, where \code{x} is a vector, usually a point from the design space.
#' It shall evaluate to the Fisher information matrix.
#' @param dDim length of \code{x}, usually the dimensionality of the design space.
#'
#' @return \code{param} returns an object of \code{class} \code{"param"}.
#' An object of class \code{"param"} is a list containing at least the following components:
#' \itemize{
#' \item fisherIf: argument
#' \item x: a row matrix of points where \code{fisherIf} has already been evaluated.
#' \item fisherI: a list of Fisher information matrices, for each row in \code{x} respectively.
#' }
#'
#' @seealso \code{\link{fisherI}}, \code{\link{update.param}}, \code{\link{Dsensitivity}}, \code{\link{getM}}, \code{\link{Defficiency}}
#'
#' @example R/examples/main.R
#'
#' @export
param = function(fisherIf, dDim) {
    r = list(fisherIf=fisherIf, x=matrix(0, nrow=0, ncol=dDim), fisherI=list())
    class(r) = 'param'
    return(r)
}


#' Update Parametric Model
#'
#' \code{update.param} evaluates the Fisher information at uncharted points and returns an updated model object.
#'
#' @param object some model.
#' @param x either a row matrix of points or a design, or a list structure of matrices or designs.
#' The number of columns/the dimensionality of the design space shall be equal to \code{ncol(object$x)}.
#' @param ... ignored.
#'
#' @return \code{update.param} returns an object of \code{class} \code{"param"}.
#' See \code{\link{param}} for its structural definition.
#'
#' @seealso \code{\link{param}}, \code{\link{design}}
#'
#' @examples ## see examples for param
#'
#' @export
update.param = function(object, x, ...) {
    mod = object # workaround for S3 requirement

    x = flatten(x)
    x = lapply(x, function(x) if (inherits(x, 'desigh')) x$x else x)
    x = do.call(rbind, x)

    if (ncol(x) != ncol(mod$x))
        stop(paste('x shall have exactly', ncol(mod$x), 'columns'))

    r = mod
    x = unique(rbind(mod$x, x)) # merge x
    idcs = seq1(nrow(mod$x) + 1, nrow(x))

    if (length(idcs) != 0) {
        xx = x[idcs,, drop=F]
        r = lapply(split(xx, seq_len(nrow(xx))), mod$fisherIf)
        fisherI = c(mod$fisherI, r)
        names(fisherI) = NULL
        mod$fisherI = fisherI
    }

    mod$x = x
    return(mod)
}


#' Build Density
#'
#' \code{buildf} builds the joint probabilty density given the marginal distributions and some copula.
#'
#' Please note that expressions are not validated.
#'
#' @param margins either \itemize{
#' \item \code{function(y, theta, ...)}, where \code{theta} is a list of parameters.
#' It shall return a column matrix of two, the probability densities and cumulative distributions.
#' \item list of pairs of expressions, named \code{"pdf"} and \code{"cdf"}, the probability density and cumulative distribution.
#' }
#' @param copula if \code{margins} is \itemize{
#' \item a function then either a copula object from package \pkg{copula} or \code{function(u, theta, ...)}, a probability density function.
#' \item a list of expressions then either a copula object from package \pkg{copula} which contains distribution expressions or an expression for the probability density.
#' }
#' @param parNames (if \code{margins} is a function and \code{copula} is a copula object) a vector of names or indices, the sequence of copula parameters in \code{theta}.
#' \code{0} or \code{""} identifies copula parameters to omit.
#'
#' @return \code{buildf} returns either \itemize{
#' \item \code{function(y, theta, ...)}, the joint probability density function, if \code{margins} is a function.
#' \item the joint probabilty density as an expression, otherwise.
#' }
#'
#' @seealso \pkg{copula}, \code{\link{expr2f}}, \code{\link{numDerivLogf}}, \code{\link{DerivLogf}}, \code{\link{fisherI}}
#'
#' @example R/examples/buildf.R
#'
#' @export
buildf = function(margins, copula, parNames=NULL) {
    tt = list(margins, copula, parNames)

    if (is.list(margins)) {
        if (inherits(copula, 'copula')) {
            att = attr(copula, 'exprdist')
            if (!is.null(att))
                copula = att$pdf
        }

        if (!is.language(copula) && !identical(copula, 1)) # special case for indepCopula
            stop('copula shall be a copula object containing distribution expressions or an expression for the probability density in this use case')

        ui = paste('u', seq1(1, length(margins)), sep='')
        uProd = parse(text=paste(ui, collapse='*'))[[1]]

        cdfs = lapply(margins, function(margin)
            substitute((a), list(a=margin$cdf)))
        names(cdfs) = ui

        pdfs = lapply(margins, function(margin)
            substitute((a), list(a=margin$pdf)))
        names(pdfs) = ui

        return(substitute((a)*b,
           list(a=substituteDirect(copula, cdfs),
                b=substituteDirect(uProd, pdfs)))
        )
    }

    if (!is.function(margins))
        stop('margins shall be a function in this use case')

    if (inherits(copula, 'copula')) {
        C = copula

        idcs = which(parNames != 0 & parNames != '')
        if (length(idcs) == 0)
            copula = function(u, theta, ...) copula::dCopula(u, C, ...)
        else {
            parNames = parNames[idcs]
            copula = function(u, theta, ...) {
                C@parameters[idcs] = as.numeric(theta[parNames])
                return(copula::dCopula(u, C, ...))
            }
        }
    }

    if (!is.function(copula))
        stop('copula shall be a copula object or a probability density function in this use case')

    return(function(y, theta, ...) {
        dp = margins(y, theta, ...)
        return(copula(dp[,2], theta) * prod(dp[,1]))
    })
}


joinLanguage = function(...) {
    x = list(...)
    x = lapply(x, function(x) if (x[[1]] == '{') as.list(x)[seq1(2, length(x))] else x)
    x = unlist(x, recursive=F, use.names=F)
    varNames = paste('x', seq1(1, length(x)), sep='')
    names(x) = varNames
    return(substituteDirect(parse(text=paste('{', paste(varNames, collapse=';'), '}', sep=''))[[1]], x))
}

withQuotes = function(x) {
    if (is.character(x))
        return(paste('\'', x, '\'', sep=''))
    return(as.character(x))
}

#' Expression To Function
#'
#' \code{expr2f} turns an expression into \code{function(y, theta, ...)}.
#'
#' @param x an expression.
#' @param map a named list of character strings defining left assignments (\code{a="b"} => \code{a <- b}).
#' @param yMap like \code{map} with \code{a=b} resolving to \code{a <- y[b]}.
#' @param thetaMap like \code{map} with \code{a=b} resolving to \code{a <- theta[[b]]}.
#'
#' @return \code{expr2f} returns \code{function(y, theta, ...)}, where \code{theta} is a list of parameters.
#' It evaluates expression \code{x}.
#'
#' @seealso \code{\link{buildf}}, \code{\link{fisherI}}
#'
#' @examples ## see examples for param
#'
#' @export
expr2f = function(x, map=NULL, yMap=NULL, thetaMap=NULL) {
    if (is.null(map))
        map = list()
    if (!is.null(yMap)) {
        yMap[] = paste('y[', lapply(yMap, withQuotes), ']', sep='')
        map = modifyList(map, yMap)
    }
    if (!is.null(thetaMap)) {
        thetaMap[] = paste('theta[[', lapply(thetaMap, withQuotes), ']]', sep='')
        map = modifyList(map, thetaMap)
    }

    map = parse(text=paste('{', paste(names(map), map, sep='=', collapse=';'), '}'))[[1]]
    return(as.function(c(alist(y=, theta=, ...=), list(joinLanguage(map, x)) )) )
}


#' Build Derivative Function for Log f
#'
#' Builds a function that evaluates to the first/second derivative of \code{log(f(y, theta, ...))} with respect to \code{theta[[i]]}/\code{theta[[i]]} and \code{theta[[j]]}.
#'
#' \pkg{numDeriv} produces \code{NaN}s if the log evaluates to (negative) \code{Inf} so you may want to specify \code{logZero} and \code{logInf}.
#'
#' @param f \code{function(y, theta, ...)}, where \code{theta} is a list of parameters.
#' A joint probability density function.
#' @param isLogf set to \code{TRUE} if \code{f} is already \code{log(f)}.
#' @param logZero the value \code{log(f)} should return if \code{f} evaluates to \code{0}.
#' @param logInf the value \code{log(f)} should return if \code{f} evaluates to \code{Inf}.
#' @param method,side,method.args see \code{\link[numDeriv]{grad}} and \code{\link[numDeriv]{hessian}} in package \pkg{numDeriv}.
#'
#' @seealso \code{\link[numDeriv]{grad}} and \code{\link[numDeriv]{hessian}} in package \pkg{numDeriv}, \code{\link{buildf}}, \code{\link{DerivLogf}}, \code{\link{fisherI}}
#'
#' @examples ## see examples for param
#'
#' @name numDerivLogf
NULL

#' @rdname numDerivLogf
#'
#' @details \code{numDerivLogf} passes \code{method}, \code{side} and \code{method.args} directly to \code{numDeriv::grad}.
#'
#' @return \code{numDerivLogf} returns \code{function(y, theta, i, ...)} which evaluates to the first derivative of \code{log(f(y, theta, ...))} with respect to \code{theta[[i]]}.
#'
#' @export
numDerivLogf = function(f, isLogf=FALSE, logZero=.Machine$double.xmin, logInf=.Machine$double.xmax/2, method='Richardson', side=NULL, method.args=list()) {
    tt = list(f, logZero, logInf, method, side, method.args)
    f = as.function(f)
    if (isLogf)
        log = function(x) x

    logf = function(theta_i, y, theta, i, ...) {
        theta[i] = theta_i
        r = log(f(y, theta, ...))
        if (is.infinite(r)) # numDeriv can't handle +-Inf
            if (r == Inf)
                return(logInf)
            else
                return(logZero)
        return(r)
    }
    r = function(y, theta, i, ...) {
        return(numDeriv::grad(logf, theta[[i]], method, side, method.args, y, theta, i, ...))
    }
    return(r)
}


#' @rdname numDerivLogf
#'
#' @details \code{numDeriv2Logf} duplicates the internals of \code{numDeriv::hessian} to gain speed.
#' The defaults for \code{method.args} are \code{list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2)}.
#'
#' @return \code{numDeriv2Logf} returns \code{function(y, theta, i, j, ...)} which evaluates to the second derivative of \code{log(f(y, theta, ...))} with respect to \code{theta[[i]]} and \code{theta[[j]]}.
#'
#' @export
numDeriv2Logf = function(f, isLogf=FALSE, logZero=.Machine$double.xmin, logInf=.Machine$double.xmax/2, method='Richardson', method.args=list()) {
    tt = list(f, logZero, logInf, method, method.args)
    f = as.function(f)
    if (isLogf)
        log = function(x) x
    method.args_ = list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2)
    method.args_[names(method.args)] = method.args

    logf = function(theta_ij, y, theta, ij, ...) {
        theta[ij] = theta_ij
        r = log(f(y, theta, ...))
        if (is.infinite(r)) # numDeriv can't handle +-Inf
            if (r == Inf)
                return(logInf)
            else
                return(logZero)
        return(r)
    }
    if (method == 'Richardson') {
        r = function(y, theta, i, j, ...) {
            if (i == j) {
                ## sole second derivative at D[2]
                return( numDeriv::genD(logf, theta[[i]], method, method.args_, y, theta, i, ...)$D[2] )
            }
            return( numDeriv::genD(logf, c(theta[[i]], theta[[j]]), method, method.args_, y, theta, c(i, j), ...)$D[4] ) # sole mixed second derivative at D[4]
        }
    } else if (method == 'complex') {
        r = function(y, theta, i, j, ...) {
            r = numDeriv::hessian(logf, c(theta[[i]], theta[[j]]), method, method.args_, y, theta, c(i, j), ...)
            if (i == j)
                return(r[1])
            return(r[2])
        }
    } else {
        stop('method not implemented.')
    }
    return(r)
}


#' Build Derivative Function for Log f
#'
#' Builds a function that evaluates to the first/second derivative of \code{log(f)} with respect to a predefined set of variables/variable combinations.
#'
#' While \code{numDerivLogf} relies on the package \pkg{numDeriv} and therefore uses finite differences to evaluate the derivatives, \pkg{DerivLogf} utilizes the package \pkg{Deriv} to build sub functions for each variable in \code{varNames}.
#' The same is true for \code{Deriv2Logf}.
#'
#' Up to version 3.6.0 of \pkg{Deriv}, \code{Deriv::Deriv} didn't recognize components or parameters accessed by \code{[}, \code{[[} or \code{$} as variables (e.g. \code{theta[["beta1"]]}).
#' Therefore it is necessary to specify mappings from \code{y} and \code{theta} to the variables in \code{f}.
#'
#' @param f an expression, a joint probability density.
#' @param varNames a character vector of variable names.
#' @param map a named list of character strings defining left assignments (\code{a="b"} => \code{a <- b}).
#' @param yMap like \code{map} with \code{a=b} resolving to \code{a <- y[b]}.
#' @param thetaMap like \code{map} with \code{a=b} resolving to \code{a <- theta[[b]]}.
#'
#' @seealso \code{\link[Deriv]{Deriv}} in package \pkg{Deriv}, \code{\link{buildf}}, \code{\link{numDerivLogf}}, \code{\link{fisherI}}
#'
#' @examples ## see examples for param
#' ## mind the gain regarding runtime compared to numDeriv
#'
#' @name DerivLogf
NULL

#' @rdname DerivLogf
#'
#' @return \code{DerivLogf} returns \code{function(y, theta, i, ...)} where \code{theta} is a list of parameters.
#' It evaluates to the first derivative of \code{log(f)} with respect to variable \code{i}.
#' Additionally the attribute \code{"d"} contains the list of sub functions.
#'
#' @export
DerivLogf = function(f, varNames, map=NULL, yMap=NULL, thetaMap=NULL) {
    logf = substitute(log((f)), list(f=f))
    d = Derivf(logf, varNames)
    d = lapply(d, expr2f, map, yMap, thetaMap)

    r = function(y, theta, i, ...) {
        return(d[[i]](y, theta, ...))
    }
    attr(r, 'd') = d
    return(r)
}

#' @rdname DerivLogf
#'
#' @return \code{Deriv2Logf} returns \code{function(y, theta, i, j, ...)} where \code{theta} is a list of parameters.
#' It evaluates to the second derivative of \code{log(f)} with respect to the variables \code{i} and \code{j}.
#' Additionally the attribute \code{"d2"} contains the list of sub functions.
#'
#' @export
Deriv2Logf = function(f, varNames, map=NULL, yMap=NULL, thetaMap=NULL) {
    logf = substitute(log((f)), list(f=f))
    d2 = Deriv2f(logf, varNames)
    d2 = lapply(d2, function(d2) lapply(d2, expr2f, map, yMap, thetaMap))

    r = function(y, theta, i, j, ...) {
        return(d2[[i]][[j]](y, theta, ...))
    }
    attr(r, 'd2') = d2
    return(r)
}


#' Fisher Information
#'
#' \code{fisherI} utilizes \code{nint_integrate} to evaluate the Fisher information.
#'
#' If \code{ff} is a list, it shall contain \code{dlogf} xor \code{d2logf}.
#'
#' @param ff either \itemize{
#' \item \code{function(y, theta, i, j, ...)} which evaluates to the inner part of the expectation integral/sum.
#' \item \code{list(f=function(y, theta, ...), d2logf=function(y, theta, i, j, ...))} (recommended)
#' \item \code{list(f=function(y, theta, ...), dlogf=function(y, theta, i, ...))}
#' }
#' where \code{f} is the joint probability density function and \code{dlogf}/\code{d2logf} is the first/second derivative of \code{log(f)} with respect to \code{theta[[i]]}/\code{theta[[i]]} and \code{theta[[j]]}.
#' @param theta the list of parameters.
#' @param parNames a vector of names or indices, the subset of parameters to use.
#' @param yspace a space, the support of \code{y}.
#' @param ... other arguments passed to \code{ff}.
#'
#' @return \code{fisherI} returns a named matrix, the Fisher information.
#'
#' @seealso \code{\link{buildf}}, \code{\link{expr2f}}, \code{\link{numDerivLogf}}, \code{\link{DerivLogf}}, \code{\link{nint_space}}, \code{\link{nint_transform}}, \code{\link{nint_integrate}}, \code{\link{param}}
#'
#' @examples ## see examples for param
#'
#' @export
fisherI = function(ff, theta, parNames, yspace, ...) {
    tt = list(ff, theta, parNames, yspace)

    i = 0
    j = 0
    if (inherits(ff, 'list')) {
        dlogf = ff[['dlogf']]
        d2logf = ff[['d2logf']]
        if ((is.null(dlogf) && is.null(d2logf)) || (!is.null(dlogf) && !is.null(d2logf)))
            stop('either dlogf xor d2logf shall be given')

        f = ff[['f']]
        if (!is.null(dlogf)) {
            g = function(y, ...) dlogf(y, theta, i, ...)*dlogf(y, theta, j, ...)*f(y, theta, ...)
            gd = function(y, ...) dlogf(y, theta, i, ...)**2 *f(y, theta, ...)
        } else {
            g = function(y, ...) -d2logf(y, theta, i, j, ...)*f(y, theta, ...)
            gd = g
        }
    } else {
        g = function(y, ...) ff(y, theta, i, j, ...)
        gd = g
    }

    n = length(parNames)
    if (n == 0)
        return(matrix(nrow=0, ncol=0))

    combs = combn(parNames, 2)
    r = matrix(0, nrow=n, ncol=n, dimnames=list(parNames, parNames))

    ## do off diagonal
    for (k in seq1(1, ncol(combs))) {
        i = combs[1, k]
        j = combs[2, k]

        r[i, j] = nint_integrate(g, yspace, ...)
        #print(j / ncol(combs))
        #print(rr)
    }

    ## do diagonal
    for (i in parNames) {
        j = i # necessary
        r[i, i] = nint_integrate(gd, yspace, ...)
        #print(j / ncol(combs))
        #print(rr)
    }

    return(mirrorMatrix(r))
}


#' Design
#'
#' \code{design} creates a custom design object.
#'
#' @param x a row matrix of points.
#' @param w a vector of weights.
#' Length shall be equal to the number of rows in \code{x} and sum shall be equal to \code{1}.
#' @param tag a list containing additional information about the design.
#'
#' @return \code{design} returns an object of \code{class} \code{"desigh"}.
#' An object of class \code{"desigh"} is a list containing at least this function's arguments.
#'
#' @seealso \code{\link{Wynn}}, \code{\link{reduce}}, \code{\link{getM}}, \code{\link{plot.desigh}}, \code{\link{Defficiency}}, \code{\link{update.param}}
#'
#' @examples ## see examples for param
#'
#' @export
design = function(x, w, tag=list()) {
    if (length(w) != nrow(x))
        stop('length of w shall be equal to the number of rows in x')
    if (!isTRUE(all.equal(sum(w), 1))) # R forced me to write this
        stop('weights shall sum to 1')

    r = list(x=x, w=w, tag=tag)
    class(r) = 'desigh'
    return(r)
}


getDAPar = function(mod, A, parNames) {
    # A = NULL || chr || int || row named mat || mat
    # parNames = NULL || chr || int

    # gather fI, fNames
    fI = NULL
    fNames = NULL

    if (length(mod$fisherI) != 0) {
        fI = mod$fisherI[[1]]

        fNames = rownames(fI)
        if (is.null(fNames)) {
            fNames = colnames(fI)
        }
    }

    # parNames
    if (is.character(parNames)) {
        if (!is.null(fNames)) {
            parNames = match(parNames, fNames)

            if (anyNA(parNames))
                stop('parNames doesn\'t correspond to the names in the first Fisher information matrix')
        }
    }
    # = NULL || chr || int

    # A
    if (is.null(A) || (is.matrix(A) && is.null(rownames(A)))) # D, raw D_A
        return(list(parNames=parNames, A=A))

    n = NULL
    aIdcs = NULL
    if (is.character(A)) {
        n = A
    } else if (is.vector(A)) {
        aIdcs = A
    } else {
        n = rownames(A)
    }

    if (!is.null(n)) { # resolve parNames
        if (is.character(parNames)) { # parNames = chr
            n2 = parNames
        } else { # parNames = NULL || int
            if (is.null(fNames))
                stop('in this case the first Fisher information matrix shall have names')

            if (is.null(parNames)) { # parNames = NULL
                n2 = fNames
            } else { # parNames = int
                n2 = fNames[parNames]
            }
        }

        aIdcs = match(n, n2)
        if (anyNA(aIdcs))
            stop('parameter names in A don\'t correspond to the subset of parameters')
    }

    if (is.null(parNames)) { # parNames = NULL
        if (is.null(fI))
            stop('in this case the model shall contain at least one Fisher information matrix')

        k = nrow(fI)
    } else { # parNames = chr || int
        k = length(parNames)
    }

    if (is.character(A) || is.vector(A)) { # D_s
        s = length(A)
        A = matrix(0, nrow=k, ncol=s)
        A[(seq1(1, s) - 1)*k + aIdcs] = 1
        return(list(parNames=parNames, A=A))
    }

    # named D_A
    nA = matrix(0, nrow=k, ncol=ncol(A))
    nA[aIdcs,] = A
    return(list(parNames=parNames, A=nA))
}

getm = function(mod, x, idcs=NULL) {
    if (nrow(x) == 0)
        return(array(dim=c(length(idcs), length(idcs), 0)))

    xidcs = rowmatch(x, mod$x)
    if (anyNA(xidcs))
        stop('model shall contain Fisher information matrices for each point. See update.param')

    fi = mod$fisherI[[xidcs[1]]]
    if (is.null(idcs) || nrow(fi) == length(idcs))
        return(simplify2array(mod$fisherI[xidcs]))

    return(vapply(xidcs, function(xidx) mod$fisherI[[xidx]][idcs, idcs], matrix(0, nrow=length(idcs), ncol=length(idcs))))
}

getM_ = function(m, w) {
    return(apply(sweep(m, 3, w, '*'), c(1, 2), sum))
}


Dsensitivity_anyNAm = 'Fisher information matrices shall not contain missing values'

#' D Sensitivity
#'
#' \code{Dsensitivity} builds a sensitivity function for the D-, D_s or D_A-optimality criterion which relies on defaults to speed up evaluation.
#' \code{Wynn} for instance requires this behaviour/protocol.
#'
#' Indices and rows of an unnamed matrix supplied to argument \code{A} correspond to the subset of parameters defined by argument \code{parNames}.
#'
#' For efficiency reasons the returned function won't complain about \emph{missing arguments} immediately, leading to strange errors.
#' Please make sure that all arguments are specified at all times.
#' This behaviour might change in future releases.
#'
#' @param A for \itemize{
#' \item D-optimality: \code{NULL}
#' \item D_s-optimality: a vector of names or indices, the subset of parameters of interest.
#' \item D_A-optimality: either \itemize{
#'   \item directly: a matrix without row names.
#'   \item indirectly: a matrix with row names corresponding to the parameters.
#'   }
#' }
#' @param parNames a vector of names or indices, the subset of parameters to use.
#' Defaults to the parameters for which the Fisher information is available.
#' @param defaults a named list of default values.
#' The value \code{NULL} is equivalent to absence.
#'
#' @return \code{Dsensitivity} returns \code{function(x=NULL, desw=NULL, desx=NULL, mod=NULL)}, the sensitivity function.
#' It's attributes contain this function's arguments.
#'
#' @seealso \code{\link{param}}, \code{\link{Wynn}}, \code{\link{plot.desigh}}
#'
#' @examples ## see examples for param
#'
#' @export
Dsensitivity = function(A=NULL, parNames=NULL, defaults=list(x=NULL, desw=NULL, desx=NULL, mod=NULL)) {
    dmod = defaults$mod
    ddesx = defaults$desx
    ddesw = defaults$desw
    dx = defaults$x

    u1 = T; u2 = T # update

    # mod, A, parNames -> parNames, A
    if (u1 && !is.null(dmod)) {
        tt = getDAPar(dmod, A, parNames)
        dnames = tt$parNames
        dA = tt$A
    } else {
        dnames = NULL
        dA = NULL
        u1 = F
    }

    # mod, parNames, desx -> desm
    if (u1 && !is.null(ddesx)) {
        ddesm = getm(dmod, ddesx, dnames)
        if (anyNA(ddesm))
            stop(Dsensitivity_anyNAm)
    } else {
        ddesm = NULL
        u1 = F
    }

    # desm, desw -> t1
    if (u1 && !is.null(ddesw)) {
        ddesM = getM_(ddesm, ddesw)
        ddesMi = solve(ddesM)
        if (is.null(dA))
            dt1 = ddesMi
        else
            dt1 = ddesMi %*% dA %*% solve(t(dA) %*% ddesMi %*% dA) %*% t(dA) %*% ddesMi
    } else {
        dt1 = NULL
        u1 = F
    }

    # mod, parNames, x -> m
    if (u1 && !is.null(dx)) {
        dm = getm(dmod, dx, dnames)
        if (anyNA(dm))
            stop(Dsensitivity_anyNAm)
    } else {
        dm = NULL
        u2 = F
    }

    r = function(x=NULL, desw=NULL, desx=NULL, mod=NULL) {
        u1 = F #; u2 = F # update

        if (is.null(mod))
            mod = dmod
        else
            u1 = T

        # mod, A, parNames -> parNames, A
        if (u1) {
            tt = getDAPar(mod, A, parNames)
            parNames = tt$parNames
            A = tt$A
        } else {
            parNames = dnames
            A = dA
        }

        if (is.null(desx))
            desx = ddesx
        else
            u1 = T

        # mod, parNames, desx -> desm
        if (u1) {
            desm = getm(mod, desx, parNames)
            if (anyNA(desm))
                stop(Dsensitivity_anyNAm)
        } else
            desm = ddesm

        if (is.null(desw))
            desw = ddesw
        else
            u1 = T

        # desm, desw -> t1
        if (u1) {
            desM = getM_(desm, desw)
            desMi = solve(desM)
            if (is.null(A))
                t1 = desMi
            else
                t1 = desMi %*% A %*% solve(t(A) %*% desMi %*% A) %*% t(A) %*% desMi
        } else
            t1 = dt1

        if (is.null(x))
            x = dx
        else
            u1 = T # u2 would be more precise

        # mod, parNames, x -> m
        if (u1) {
            m = getm(mod, x, parNames)
            if (anyNA(m))
                stop(Dsensitivity_anyNAm)
        } else
            m = dm

        return(apply(m, 3, function(m, t1) sum(diag(t1 %*% m)), t1))
    }
    attributes(r) = list(A=A, parNames=parNames, defaults=defaults)

    return(r)
}


#' Wynn
#'
#' \code{Wynn} finds an optimal design using some sensitivity function and a Wynn-algorithm.
#'
#' See \code{\link{Dsensitivity}} and it's return value for a reference implementation of a function complying with the requirements for \code{sensF}.
#'
#' The algorithm starts from a uniform weight design.
#' In each iteration weight is redistributed to the point which has the highest sensitivity.
#' Sequence: \code{1/i}.
#' The algorithm stops when all sensitivities are below a specified tolerance level or the maximum number of iterations is reached.
#'
#' @param sensF \code{function(x=NULL, desw=NULL, desx=NULL, mod=NULL)}, a sensitivity function.
#' It's attribute \code{"defaults"} shall contain identical \code{x} and \code{desx}, and \code{sensF(desw=w)} shall return sensitivities corresponding to each point in \code{x}.
#' @param tol the tolerance level regarding the sensitivities.
#' @param maxIter the maximum number of iterations.
#'
#' @return \code{Wynn} returns an object of \code{class} \code{"desigh"}.
#' See \code{\link{design}} for its structural definition.
#'
#' @references Wynn, Henry P. (1970) The Sequential Generation of D-Optimum Experimental Designs.
#' \emph{The Annals of Mathematical Statistics}, 41(5):1655-1664.
#'
#' @seealso \code{\link{Dsensitivity}}, \code{\link{design}}
#'
#' @examples ## see examples for param
#'
#' @export
Wynn = function(sensF, tol, maxIter=1e4) {
    tag = list(Wynn=list(sensF=sensF, tol=tol, maxIter=maxIter))

    defaults = attr(sensF, 'defaults')
    x = defaults$x
    if (!isTRUE(all.equal(x, defaults$desx)))
        stop('sensitivity defaults for x and desx shall be equal')

    if (nrow(x) == 0)
        return(design(x, numeric(0), tag=tag))

    n = nrow(x)
    w = rep(1/n, n)
    tolBreak = F

    for (iIter in seq1(1, maxIter)) {
        sens = sensF(desw=w)

        maxIdx = which(sens == max(sens))
        if (length(maxIdx) != 1)
            maxIdx = sample(maxIdx, 1)

        dw = 1 / (iIter + 1)
        w = w * (1 - dw)
        w[maxIdx] = 0
        w[maxIdx] = 1 - sum(w) # equal to 'w[maxIdx] + dw'

        mSens = sens[maxIdx]
        if (mSens <= tol) {
            tolBreak = T
            break
        }
    }

    names(sens) = NULL
    tag$Wynn$tolBreak = tolBreak
    tag$Wynn$nIter = iIter
    return(design(x, w, tag=tag))
}


wPoint = function(x, w) {
    ## x = row matrix
    ## w = vector
    ## nrow(x) == length(w)
    return( apply(sweep(x, 1, w, '*'), 2, sum) / sum(w) )
}

#' Reduce Design
#'
#' \code{reduce} drops insignificant points and merges points in a certain neighbourhood.
#'
#' @param des some design.
#' @param distMax maximum euclidean distance between points to be merged.
#' @param wMin minimum weight a point shall have to be considered significant.
#'
#' @return \code{reduce} returns an object of \code{class} \code{"desigh"}.
#' See \code{\link{design}} for its structural definition.
#'
#' @seealso \code{\link{design}}
#'
#' @examples ## see examples for param
#'
#' @export
reduce = function(des, distMax, wMin=1e-6) {
    x = des$x
    w = des$w

    m = wMin <= w
    x = x[m,, drop=F]
    w = w[m]
    cl = clusterPeak(x, w, distMax)

    wg = split(w, cl)
    rx = do.call(rbind, mapply(function(x, w) wPoint(x, w), split(as.data.frame(x), cl), wg, SIMPLIFY=F))
    if (is.null(rx))
        rx = matrix(nrow=0, ncol=ncol(x))
    dimnames(rx) = NULL
    rw = vapply(wg, sum, 0)
    rw = rw / sum(rw)
    names(rw) = NULL

    ord = roworder(rx)
    rx = rx[ord,, drop=F]
    rw = rw[ord]

    tag = list(reduce=list(des=des, distMax=distMax, wMin=wMin))
    return(design(rx, rw, tag=tag))
}


#' Get Fisher Information
#'
#' \code{getM} returns the Fisher information corresponding to some model and some design.
#'
#' @param mod some model.
#' @param des some design.
#'
#' @return \code{getM} returns a named matrix, the Fisher information.
#'
#' @seealso \code{\link{param}}, \code{\link{design}}
#'
#' @examples ## see examples for param
#'
#' @export
getM = function(mod, des) {
    m = getm(mod, des$x)
    return(getM_(m, des$w))
}


## from http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
## Adds transparency to colours
add.alpha <- function(col, alpha=1){
    if (missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))
}

getBaseDesign = function(des) {
    tag = des$tag
    if (!is.null(tag$Wynn))
        return(des)
    if (!is.null(tag$reduce))
        return(getBaseDesign(tag$reduce$des))
    return(NULL)
}

#' Plot Design
#'
#' \code{plot.desigh} creates a one-dimensional design plot, optionally together with a specified sensitivity curve.
#' If the design space has additional dimensions, the design is projected on a specified margin.
#'
#' @param x some design.
#' @param sensx (optional) a row matrix of points.
#' @param sens (optional) either a vector of sensitivities or a sensitivity function.
#' The latter shall rely on defaults, see \code{\link{Dsensitivity}} for details.
#' @param sensTol (optional) a single numeric.
#' Adds a horizontal line at this sensitivity level.
#' @param ... other arguments passed to plot.
#' @param margins a vector of indices, the dimensions to project on.
#' Defaults to \code{1}.
#' @param desSens if \code{TRUE} and \code{sens} is not specified then the sensitivity function which potentially was used in \code{Wynn} is taken as \code{sens}.
#' @param sensPch either a character vector of point 'characters' to add to the sensitivity curve or \code{NULL}.
#' @param sensArgs a list of arguments passed to draw calls related to the sensitivity.
#'
#' @references uses add.alpha from \url{http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html}
#'
#' @seealso \code{\link{design}}, \code{\link{Dsensitivity}}
#'
#' @examples ## see examples for param
#'
#' @export
plot.desigh = function(x, sensx=NULL, sens=NULL, sensTol=NULL, ..., margins=NULL, desSens=T, sensPch='+', sensArgs=list()) {
    points = function(..., axes) graphics::points(...)
    abline = function(..., axes) graphics::abline(...)
    axis = function(..., axes) graphics::axis(...)
    mtext = function(..., axes) graphics::mtext(...)

    des = x # workaround for S3 requirement

    args = list(...)

    if (is.null(margins))
        margins = 1
        #margins = 1:(ncol(des$model$x))
    if (1 < length(margins))
        stop('not yet implemented')

    ## marginal projection
    x = des$x
    w = des$w
    idcs = split(seq1(1, nrow(x)), lapply(margins, function(margin) x[, margin]), drop=T)
    x = x[sapply(idcs, function(i) i[1]), margins, drop=F]
    w = sapply(idcs, function(idcs) sum(w[idcs]))
    ord = roworder(x)
    x = x[ord,, drop=F]
    w = w[ord]

    ## lookup design sensF
    if (is.null(sens) && desSens) {
        d = getBaseDesign(des)
        sens = d$tag$Wynn$sensF
        if (is.null(sensx) && !is.null(sens)) {
            sensx = attr(sens, 'defaults')$x
            sens = sens(desw=d$w)
        }
    }

    ## prepare sensitivity
    if (!is.null(sens)) {
        if (is.function(sens)) {
            if (is.null(sensx)) {
                sensx = attr(sens, 'defaults')$x
                sens = sens(desw=des$w, desx=des$x)
            } else
                sens = sens(x=sensx, desw=des$w, desx=des$x)
        }

        if (is.null(sensx))
            stop('if sens is a vector, sensx shall be specified')

        ## marginal projection
        idcs = split(seq1(1, nrow(sensx)), lapply(margins, function(margin) sensx[, margin]), drop=T)

        sensx = sensx[sapply(idcs, function(i) i[1]), margins, drop=F]
        sens = sapply(idcs, function(idcs) max(sens[idcs]))

        ord = roworder(sensx)
        sensx = sensx[ord,, drop=F]
        sens = sens[ord]

        ## add axis margin
        mar = par('mar')
        mar[4] = mar[2]
        par(mar=mar)
    } else {
        ## remove axis margin
        mar = par('mar')
        mar[4] = 2 + 0.1 # Warning, hard coded
        par(mar=mar)
    }

    if (length(margins) == 1) {
        xlab = colnames(x)
        if (is.null(xlab))
            xlab = paste('x[, c(', toString(margins), ')]', sep='')

        margs = modifyList(list(NA, xlim=range(x), ylim=c(0, 1), xlab=xlab, ylab='weight'), args)
        do.call(plot, margs)

        xlim = margs$xlim; ylim = margs$ylim

        if (!is.null(sens)) {
            par(new=T)

            dylim = range(sens)
            if (0 < dylim[1])
                dylim = c(0, dylim[2])
            else if (dylim[2] < 0)
                dylim = c(dylim[1], 0)
            if (isTRUE(sensTol < dylim[1]))
                dylim = c(sensTol, dylim[2])
            else if (isTRUE(dylim[2] < sensTol))
                dylim = c(dylim[1], sensTol)

            defaultYlab = ifelse(ncol(des$x) == 1, 'sensitivity', 'maximum sensitivity')
            margs = modifyList(list(sensx, sens, type='l', ylim=dylim, ylab=defaultYlab), sensArgs)
            ylab = margs$ylab
            margs = modifyList(margs, list(xlim=xlim, axes=F, xlab='', ylab=''))
            do.call(plot, margs)

            if (!is.null(sensTol)) {
                margs = modifyList(list(h=sensTol, col='black', lty=2), sensArgs)
                #margs$col = add.alpha(margs$col, 0.33)
                do.call(abline, margs)
            }

            if (!isTRUE(sensArgs$axes == F)) {
                margs = modifyList(list(4), sensArgs)
                if (!is.null(margs$col) && is.null(margs$col.axis))
                    margs$col.axis = margs$col
                do.call(axis, margs)
            }

            margs = modifyList(list(ylab, side=4, line=3), sensArgs)
            do.call(mtext, margs)

            if (!(is.null(sensPch) || identical(sensPch, ''))) {
                alpha = des$w / max(des$w)
                idcs = which(1/256 < alpha)
                px = des$x[idcs,]
                p = approx(sensx, sens, px)

                margs = modifyList(list(p$x, p$y, pch=sensPch, col='black'), args)
                margs$col = add.alpha(margs$col, alpha[idcs])
                do.call(points, margs)
            }
        }

        par(new=T)
        margs = modifyList(list(x, w, xlim=xlim, ylim=ylim, type='h'), args)
        margs = modifyList(margs, list(axes=F, xlab='', ylab=''))
        do.call(plot, margs)
    }
}


#' D Efficiency
#'
#' \code{Defficiency} computes the D-, D_s or D_A-efficiency measure for some design with respect to some reference design.
#'
#' Indices supplied to argument \code{A} correspond to the subset of parameters defined by argument \code{parNames}.
#'
#' D efficiency is defined as
#' \deqn{\left(\frac{\left|M(\xi,\bar{\theta})\right|}{\left|M(\xi^{*},\bar{\theta})\right|}\right)^{1/n}}{( det(M(\xi, \theta)) / det(M(\xi*, \theta)) )**(1/n)}
#' and D_A efficiency as
#' \deqn{\left(\frac{\left|A^{T}M(\xi^{*},\bar{\boldsymbol{\theta}})^{-1}A\right|}{\left|A^{T}M(\xi,\bar{\boldsymbol{\theta}})^{-1}A\right|}\right)^{1/s}}{( det(t(A) \%*\% solve(M(\xi*, \theta)) \%*\% A) / det(t(A) \%*\% solve(M(\xi, \theta)) \%*\% A) )**(1/s)}
#'
#' @param des some design.
#' @param ref some design, the reference.
#' @param mod some model.
#' @param A for \itemize{
#' \item D-efficiency: \code{NULL}
#' \item D_s-efficiency: a vector of names or indices, the subset of parameters of interest.
#' \item D_A-efficiency: either \itemize{
#'   \item directly: a matrix without row names.
#'   \item indirectly: a matrix with row names corresponding to the parameters.
#'   }
#' }
#' @param parNames a vector of names or indices, the subset of parameters to use.
#' Defaults to the parameters for which the Fisher information is available.
#'
#' @return \code{Defficiency} returns a single numeric.
#'
#' @seealso \code{\link{design}}, \code{\link{param}}
#'
#' @examples ## see examples for param
#'
#' @export
Defficiency = function(des, ref, mod, A=NULL, parNames=NULL) {
    tt = getDAPar(mod, A, parNames)
    parNames = tt$parNames
    A = tt$A

    m = getm(mod, des$x, parNames)
    M = getM_(m, des$w)

    m = getm(mod, ref$x, parNames)
    Mref = getM_(m, ref$w)

    if (is.null(A))
        return( (det(M) / det(Mref))**(1/nrow(M)) )

    return( (det(t(A) %*% solve(Mref) %*% A) / det(t(A) %*% solve(M) %*% A) )**(1/ncol(A)) )
}

