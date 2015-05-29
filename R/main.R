

#' Create Parametric Model
#'
#' \code{parm} creates an initial model object.
#' Unlike other model statements this function does not perform any computation.
#'
#' @param fisherIf \code{function(x, ...)}, where \code{x} is a numeric vector, usually a point from the design space. It shall evaluate to the Fisher information.
#' @param dDim length of \code{x}, usually the dimensionality of the design space.
#'
#' @return \code{parm} returns an object of \code{class} \code{"parm"}.
#' An object of class \code{"parm"} is a list containing the following components:
#'
#' fisherIf: argument
#'
#' dDim: argument
#'
#' x: row matrix of points where \code{fisherIf} has already been evaluated.
#'
#' fisherI: list of Fisher information matrices for each row of \code{x} respectively.
#'
#' @seealso \code{\link{update.parm}}
#'
#' @export
parm = function(fisherIf, dDim) {
    r = list(fisherIf=fisherIf, dDim=dDim, x=matrix(nrow=0, ncol=dDim), fisherI=list())
    class(r) = 'parm'
    return(r)
}


#' Update Parametric Model
#'
#' \code{update.parm} evaluates the Fisher information if necessary and returns the updated model object.
#'
#' @param mod an object of class \code{"parm"}
#' @param x a row matrix of points to evaluate the Fisher information at.
#' The number of columns shall be equal to \code{mod$dDim}.
#' @return \code{update.parm} returns an object of \code{class} \code{"parm"}.
#'
#' @seealso \code{\link{parm}}
#'
#' @export
update.parm = function(mod, x) {
    if (ncol(x) != mod$dDim)
        stop(paste('x shall have exactly', mod$dDim, 'columns'))

    x = unique(rbind(mod$x, x)) # merge x
    idcs = seqi(nrow(mod$x) + 1, nrow(x))

    if (length(idcs) != 0) {
        xx = x[idcs,,drop=F]
        r = lapply(split(xx, row(xx)), mod$fisherIf)
        fisherI = c(mod$fisherI, r)
        names(fisherI) = NULL
    }

    mod$x = x
    mod$fisherI = fisherI

    return(mod)
}


# from http://adv-r.had.co.nz/Computing-on-the-language.html#substitute
substitute_q <- function(x, env) {
    call <- substitute(substitute(y, env), list(y = x))
    eval(call)
}


#' Build Density
#'
#' \code{buildf} builds an object that evaluates to the density of a random vector given the marginal distributions and some copula.
#'
#' @param margins either \itemize{
#' \item \code{function(y, theta, ...)}, where \code{theta} is a list of parameters. It shall return a matrix with densities in the first column and cumulative densities in the second.
#' \item list of expressions each of which contains the PDF and CDF accessible by \code{$pdf} and \code{$cdf}. See examples below.
#' }
#' @param copula a copula object from package \pkg{copula}.
#' @param alphaIdcs an integer vector specifying the positions of copula parameters in \code{theta}.
#' Irrelevant if \code{margins} is a list.
#'
#' @return \code{buildf} returns either \itemize{
#' \item \code{function(y, theta, ...)}, the joint PDF, if \code{margins} is a function.
#' \item the joint PDF as an expression built upon the copula PDF, the marginal PDFs and CDFs, otherwise.
#' }
#'
#' @details If \code{buildf} should build an expression, the copula needs to contain the probability density as an expression.
#' Note also that there is no check whether or not the final expression is valid.
#'
#' @seealso \pkg{copula}, \code{\link{expr2f}}, \code{\link{fisherI}}, \code{\link{numDerivLogf}}, \code{\link{DerivLogf}}
#'
#' @example examples/buildf.R
#'
#' @export
buildf = function(margins, copula, alphaIdcs=NULL) {
    tt = list(margins, copula, alphaIdcs)
    if (is.function(margins)) {
        if (length(alphaIdcs) == 0) {
            r = function(y, theta, ...) {
                dp = margins(y, theta, ...)
                return(copula::dCopula(dp[,2], copula) * prod(dp[,1]))
            }
        } else {
            r = function(y, theta, ...) {
                dp = margins(y, theta, ...)
                copula@parameters = as.numeric(theta[alphaIdcs])
                return(copula::dCopula(dp[,2], copula) * prod(dp[,1]))
            }
        }
    } else { # list of expressions
        exprdist = attr(copula, 'exprdist')
        if (is.null(exprdist))
            stop('copula does not have closed form distribution expressions')

        margins = lapply(margins, function(margin) lapply(margin, function(e)  substitute((e), list(e=e)) )) # wrap expressions in ()
        names_ = paste('u', 1:length(margins), sep='')
        p = lapply(margins, function(margin) margin$cdf)
        names(p) = names_
        d = lapply(margins, function(margin) margin$pdf)
        names(d) = names_

        r = substitute((a)*b, list(a=substitute_q(exprdist$pdf, p), b=substitute_q(parse(text=paste(names_, collapse='*'))[[1]], d)))
    }
    return(r)
}


mergeLanguage = function(...) {
    x = list(...)
    x = lapply(x, function(x) if (x[[1]] == '{') as.list(x)[seqi(2, length(x))] else x)
    x = unlist(x, recursive=F, use.names=F)
    names_ = paste('x', seqi(1, length(x)), sep='')
    names(x) = names_
    return(substitute_q(parse(text=paste('{', paste(names_, collapse=';'), '}', sep=''))[[1]], x))
}

withQuotes = function(x) {
    if (is.character(x))
        return(paste('\'', x, '\'', sep=''))
    return(as.character(x))
}

#' Expression To Function
#'
#' Turns an expression into \code{function(y, theta, ...)}.
#'
#' @param x an expression.
#' @param map a named list of names defining left assignments (\code{a="b"} := \code{a <- b}).
#' @param yMap like \code{map} but items \code{a=1} resolve to \code{a <- y[1]}.
#' @param thetaMap like \code{yMap} for \code{theta} with single element access \code{[[}.
#'
#' @return \code{expr2f} returns \code{function(y, theta, ...)}, where \code{theta} is a list of parameters. It evaluates expression \code{x}.
#'
#' @seealso \code{\link{buildf}}, \code{\link{fisherI}}
#'
#' @example examples/expr2f.R
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
    return(as.function(c(alist(y=, theta=, ...=), list(mergeLanguage(map, x)) )) )
}


#' Build Derivative Function for Log f
#'
#' Builds a function that evaluates to the first or second derivative of \code{log(f(y, theta, ...))} with respect to \code{theta[[i]]} or \code{theta[[i]]} and \code{theta[[j]]}.
#'
#' \pkg{numDeriv} produces \code{NaN}s if the log evaluates to (negative) \code{Inf} so you may want to specify \code{logZero} and \code{logInf}.
#'
#' @param f \code{function(y, theta, ...)}, where \code{theta} is a list of parameters.
#' A joint PDF.
#' @param logZero the value \code{log(f)} should return if \code{f} evaluates to 0. See details below.
#' @param logInf the value \code{log(f)} should return if \code{f} evaluates to \code{Inf}.
#' @param method see \pkg{numDeriv}
#' @param method.args see \pkg{numDeriv}
#'
#' @seealso \code{\link{buildf}}, \code{\link{DerivLogf}}, \code{\link{fisherI}}, \code{\link[numDeriv]{grad}} and \code{\link[numDeriv]{hessian}} in package \pkg{numDeriv}
#'
#' @example examples/numDerivLogf.R
#'
#' @name numDerivLogf
NULL

#' @rdname numDerivLogf
#'
#' @return \code{numDerivLogf} returns \code{function(y, theta, i, ...)} which evaluates to the first derivative of \code{log(f(y, theta, ...))} with respect to \code{theta[[i]]}.
#'
#' @export
numDerivLogf = function(f, logZero=.Machine$double.xmin, logInf=.Machine$double.xmax/2, method='Richardson', method.args=list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=F)) {
    tt = list(f, logZero, logInf, method, method.args)
    f = as.function(f)

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
        return(numDeriv::grad(logf, theta[[i]], method, NULL, method.args, y, theta, i, ...))
    }
    return(r)
}


#' @rdname numDerivLogf
#'
#' @return \code{numDeriv2Logf} returns \code{function(y, theta, i, j, ...)} which evaluates to the second derivative of \code{log(f(y, theta, ...))} with respect to \code{theta[[i]]} and \code{theta[[j]]}.
#'
#' @export
numDeriv2Logf = function(f, logZero=.Machine$double.xmin, logInf=.Machine$double.xmax/2, method='Richardson', method.args=list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=F)) {
    tt = list(f, logZero, logInf, method, method.args)
    f = as.function(f)

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
                return( numDeriv::genD(logf, theta[[i]], method, method.args, y, theta, i, ...)$D[2] )
                # sole second derivative at D[2]
            }
            return( numDeriv::genD(logf, c(theta[[i]], theta[[j]]), method, method.args, y, theta, c(i, j), ...)$D[4] ) # sole mixed second derivative at D[4]
        }
    } else if (method == 'complex') {
        r = function(y, theta, i, j, ...) {
            r = numDeriv::hessian(logf, c(theta[[i]], theta[[j]]), method, method.args, y, theta, c(i, j), ...)
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
#' Builds a function that evaluates to the first or second derivative of \code{log(f)} with respect to one or any two variables specified.
#'
#' While \code{numDerivLogf} relies on \pkg{numDeriv} and therefore uses finite differences to evaluate the derivatives, \code{DerivLogf} utilizes \code{Deriv} to build sub functions for each variable in argument \code{names}.
#' The same is true for \code{Deriv2Logf}.
#'
#' \code{Deriv} won't recognize parameters accessed by \code{[[} (e.g. \code{theta[["beta1"]]}).
#' Therefore you need to specify mappings from \code{y} and \code{theta} to the variables used in \code{f}.
#' See arguments above and examples below.
#'
#' @param f a joint PDF as an expression.
#' @param names a character vector of variable names.
#' @param map a named list of names defining left assignments (\code{a="b"} := \code{a <- b}). See details below.
#' @param yMap like \code{map} but items \code{a=1} resolve to \code{a <- y[1]}.
#' @param thetaMap like \code{yMap} for \code{theta} with single element access \code{[[}.
#'
#' @seealso \code{\link{buildf}}, \code{\link{numDerivLogf}}, \code{\link{fisherI}}, \code{\link[Deriv]{Deriv}} in package \pkg{Deriv}
#'
#' @example examples/DerivLogf.R
#'
#' @name DerivLogf
NULL

#' @rdname DerivLogf
#'
#' @return \code{DerivLogf} returns \code{function(y, theta, i, ...)} where \code{theta} is a list of parameters.
#' It evaluates to the first derivative of \code{log(f)} with respect to variable \code{i}.
#' Additionally the attribute \code{"d"} contains the list of derivative functions.
#'
#' @export
DerivLogf = function(f, names, map=NULL, yMap=NULL, thetaMap=NULL) {
    logf = substitute(log((f)), list(f=f))
    d = Derivf(logf, names)
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
#' Additionally the attribute \code{"d2"} contains the list of all pairwise derivative functions.
#'
#' @export
Deriv2Logf = function(f, names, map=NULL, yMap=NULL, thetaMap=NULL) {
    logf = substitute(log((f)), list(f=f))
    d2 = Deriv2f(logf, names)
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
#' Either \code{dlogf} xor \code{d2logf} shall be given.
#'
#' @param f \code{function(y, theta, ...)}, a joint PDF.
#' @param theta a list of parameters.
#' @param names a vector of names or indices, the subset of parameters.
#' @param dlogf \code{function(y, theta, i, ...)} which evaluates to the first derivative of \code{log(f)} with respect to \code{theta[[i]]}.
#' @param d2logf \code{function(y, theta, i, j, ...)} which evaluates to the second derivative of \code{log(f)} with respect to \code{theta[[i]]} and \code{theta[[j]]}.
#' @param yspace the support of \code{y}. See \code{\link{nint_space}}.
#' @param ... other arguments passed to \code{d2logf}.
#' @param yispaces see \code{\link{nint_integrate}}.
#' @param transformInf see \code{\link{nint_integrate}}.
#'
#' @return \code{fisherI} returns a named matrix, the Fisher information.
#'
#' @seealso \code{\link{buildf}}, \code{\link{numDerivLogf}}, \code{\link{DerivLogf}}, \code{\link{nint_space}}, \code{\link{nint_integrate}}
#'
#' @examples TODO
#'
#' @export
fisherI = function(f, theta, names, dlogf=NULL, d2logf=NULL, yspace=NULL, ..., yispaces=NULL, transformInf=F) {
    tt = list(f, theta, names, dlogf, d2logf, yspace, yispaces, transformInf)

    n = length(names)
    if (n == 0)
        return(matrix(nrow=0, ncol=0))

    if ((is.null(dlogf) && is.null(d2logf)) || (!is.null(dlogf) && !is.null(d2logf)))
        stop('either dlogf xor d2logf shall be given')
    first = !is.null(dlogf)

    if (!is.null(yspace))
        yispaces = c(nint_ispaces(yspace), yispaces)

    transformInf_ = transformInf

    combs = combn(names, 2)
    r = matrix(0, nrow=n, ncol=n, dimnames=list(names, names))

    # prepare off diagonal
    if (first) {
        g = function(y, ...) dlogf(y, theta, i, ...)*dlogf(y, theta, j, ...)*f(y, theta, ...)
    } else {
        g = function(y, ...) d2logf(y, theta, i, j, ...)*f(y, theta, ...)
    }

    # do off diagonal
    for (k in 1:ncol(combs)) {
        i = combs[1, k]
        j = combs[2, k]
        
        r[i, j] = nint_integrate(g, NULL, ..., ispaces=yispaces, transformInf=transformInf_)
        #print(j / ncol(combs))
        #print(rr)
    }

    # prepare diagonal
    if (first) {
        g = function(y, ...) dlogf(y, theta, i, ...)**2 *f(y, theta, ...)
    } else {
        g = function(y, ...) d2logf(y, theta, i, i, ...)*f(y, theta, ...)
    }

    # do diagonal
    for (i in 1:n) {
        r[i, i] = nint_integrate(g, NULL, ..., ispaces=yispaces, transformInf=transformInf_)
        #print(j / ncol(combs))
        #print(rr)
    }

    r = mirrorMatrix(r)
    if (!first)
        r = -r
    return(r)
}
