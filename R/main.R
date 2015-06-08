

#' Create Parametric Model
#'
#' \code{parm} creates an initial model object.
#' Unlike other model statements this function does not perform any computation.
#'
#' @param fisherIf \code{function(x, ...)}, where \code{x} is a numeric vector, usually a point from the design space.
#' It shall evaluate to the Fisher information.
#' @param dDim length of \code{x}, usually the dimensionality of the design space.
#'
#' @return \code{parm} returns an object of \code{class} \code{"parm"}.
#' An object of class \code{"parm"} is a list containing the following components:
#' \itemize{
#' \item fisherIf: argument
#' \item dDim: argument
#' \item x: row matrix of points where \code{fisherIf} has already been evaluated.
#' \item fisherI: list of Fisher information matrices for each row of \code{x} respectively.
#' }
#'
#' @seealso \code{\link{update.parm}}, \code{\link{FedorovWynn}}
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
#' @param mod an object of class \code{"parm"}.
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

    r = mod
    x = unique(rbind(mod$x, x)) # merge x
    idcs = seqi(nrow(mod$x) + 1, nrow(x))

    if (length(idcs) != 0) {
        xx = x[idcs,, drop=F]
        r = lapply(split(xx, row(xx)), mod$fisherIf)
        fisherI = c(mod$fisherI, r)
        names(fisherI) = NULL
        mod$fisherI = fisherI
    }

    mod$x = x
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
#' \item \code{function(y, theta, ...)}, where \code{theta} is a list of parameters.
#' It shall return a matrix with densities in the first column and cumulative densities in the second.
#' \item list of expressions each of which contains the PDF and CDF accessible by \code{$pdf} and \code{$cdf}.
#' See examples below.
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
        names_ = paste('u', seqi(1, length(margins)), sep='')
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
#' \code{expr2f} turns an expression into \code{function(y, theta, ...)}.
#'
#' @param x an expression.
#' @param map a named list of names defining left assignments (\code{a="b"} := \code{a <- b}).
#' @param yMap like \code{map} but items \code{a=1} resolve to \code{a <- y[1]}.
#' @param thetaMap like \code{yMap} for \code{theta} with single element access \code{[[}.
#'
#' @return \code{expr2f} returns \code{function(y, theta, ...)}, where \code{theta} is a list of parameters.
#' It evaluates expression \code{x}.
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
#' @param logZero the value \code{log(f)} should return if \code{f} evaluates to 0.
#' See details below.
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
#' @param map a named list of names defining left assignments (\code{a="b"} := \code{a <- b}).
#' See details below.
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
#' @param yspace the support of \code{y}.
#' See \code{\link{nint_space}}.
#' @param ... other arguments passed to \code{d2logf}.
#' @param yispaces see \code{\link{nint_integrate}}.
#' @param transformInf see \code{\link{nint_integrate}}.
#'
#' @return \code{fisherI} returns a named matrix, the Fisher information.
#'
#' @seealso \code{\link{buildf}}, \code{\link{numDerivLogf}}, \code{\link{DerivLogf}}, \code{\link{nint_space}}, \code{\link{nint_integrate}}
#'
#' @examples #TODO
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
    for (k in seqi(1, ncol(combs))) {
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



design = function(mod, x, w, sens, args=list(), adds=list()) {
    r = list(model=mod, x=x, w=w, sens=sens, args=args, adds=adds)
    class(r) = 'design'
    return(r)
}

getM = function(m, w) {
    return(apply(sweep(m, 3, w, '*'), c(1, 2), sum))
}

sensD = function(m, Mi, ...) {
    return(apply(m, 3, function(m, Mi) sum(diag(Mi %*% m)), Mi))
}

#checkDsA = function(A) {
    #if (ncol(A) < nrow(A))
        #stop('\'A\' shall at least have the same number of columns as rows')
    #if (!identical(A[, seqi(1, nrow(A)), drop=F], diag(nrow(A)) ))
        #stop('the left part of \'A\' shall be the identity matrix')
    #if (any(A[, seqi(nrow(A) + 1, ncol(A))] != 0))
        #stop('the right part of \'A\' shall be zero')
#}

sensDs = function(m, Mi, A, ...) {
    t1 = Mi %*% A %*% solve(t(A) %*% Mi %*% A) %*% t(A) %*% Mi
    return(apply(m, 3, function(m, t1) sum(diag(t1 %*% m)), t1))
}

getIdcs = function(names, mod) {
    if (is.numeric(names))
        return(unique(names))

    if (length(mod$fisherI) == 0)
        stop('model shall contain at least one Fisher information matrix')
    tt = mod$fisherI[[1]]

    if (is.null(names))
        return( seqi(1, nrow(tt)) )

    nn = rownames(tt)
    if (is.null(nn)) {
        nn = colnames(tt)
        if (is.null(nn)) {
            stop('first Fisher information matrix shall contain row or column names')
        }
    }

    return( unique(match(names, nn)) )
}

getA = function(idcs, n) {
    s = length(idcs)
    r = matrix(0, nrow=n, ncol=s)
    r[(seqi(1, s) - 1)*n + idcs] = 1
    return(r)
}

#' Fedorov Wynn Algorithm
#'
#' \code{FedorovWynn} computes a D- or Ds-optimal design using the Fedorov-Wynn algorithm.
#'
#' TODO shortly describe what is done and why (reference to paper).
#'
#' @param mod an object of \code{class} \code{"parm"}.
#' @param names a vector of names or indices, the subset of parameters to optimize for.
#' @param tolAbs the absolute tolerance for the upper bound of the sensitivity.
#' @param tolRel the relative tolerance with respect to the number of parameters.
#' @param maxIter the maximum number of iterations.
#'
#' @return \code{FedorovWynn} returns an object of \code{class} \code{"design"}.
#' An object of class \code{"design"} is a list containing the following components:
#' \itemize{
#' \item mod: argument
#' \item x: row matrix of design points, in this case \code{mod$x}.
#' \item w: a numeric vector of weights for each design point.
#' \item sens: a numeric vector of sensitivity for each design point.
#' \item args: a list of arguments.
#' \item adds: a list of additional (runtime) information.
#' }
#'
#' @examples #TODO
#'
#' @seealso \code{\link{parm}}, \code{\link{reduce}}
#'
#' @export
FedorovWynn = function(mod, names=NULL, tolAbs=Inf, tolRel=1e-4, maxIter=1e4) {
    args = list(FedorovWynn=list(names=names, tolAbs=tolAbs, tolRel=tolRel, maxIter=maxIter))

    if (nrow(mod$x) == 0)
        return(design(mod, mod$x, numeric(0), numeric(0), args=args))

    m = simplify2array(mod$fisherI)

    n = dim(m)[1]
    idcs = getIdcs(names, mod)
    if ( identical(idcs, seqi(1, n)) ) {
        sensF = sensD
        A = NULL
        target = n
    } else {
        sensF = sensDs
        s = length(idcs)
        A = getA(idcs, n)
        target = s
    }

    tolAbs_ = min(tolAbs, tolRel * target)
    n = dim(m)[3]
    w = rep(1/n, n)

    for (iIter in seqi(1, maxIter)) {
        M = getM(m, w)
        Mi = solve(M) # TODO might break

        sens = sensF(m, Mi, A)

        maxIdx = which.max(sens)
        d = sens[maxIdx] - target
        if (d < tolAbs_)
            break

        dw = 1 / (iIter + 1)
        w = w * (1 - dw)
        w[maxIdx] = 0
        w[maxIdx] = 1 - sum(w) # equal to 'w[maxIdx] + dw'
    }

    base::names(sens) = NULL
    adds = list(FedorovWynn=list(tolBreak=d < tolAbs_, nIter=iIter))
    return(design(mod, mod$x, w, sens, args=args, adds=adds))
}


wPoint = function(x, w) {
    # x = row matrix
    # w = vector
    # nrow(x) == length(w)
    return( apply(sweep(x, 1, w, '*'), 2, sum) / sum(w) )
}


#' Reduce Design
#'
#' \code{reduce} drops insignificant design points and merges design points in a certain neighbourhood.
#'
#' @param des an object of \code{class} \code{"design"}.
#' @param distMax maximum euclidean distance between points to be merged.
#' @param wMin minimum weight a significant design point shall have.
#'
#' @return \code{reduce} returns an object of \code{class} \code{"design"}.
#' See \code{FedorovWynn} for its structural definition.
#' The sensitivity is \code{NA} for every design point.
#'
#' @seealso \code{\link{FedorovWynn}}, \code{\link{update.design}}
#'
#' @examples #TODO
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

    ord = orderMatrix(rx)
    rx = rx[ord,, drop=F]
    rw = rw[ord]
    rargs = des$args
    rargs$reduce = list(des=des, distMax=distMax, wMin=wMin)
    return(design(des$model, rx, rw, rep(NA, nrow(rx)), rargs, des$adds))
}


getm = function(des) {
    idcs = indexMatrix(des$model$x, des$x)
    if (anyNA(idcs))
        stop('model shall contain Fisher information matrices for each point in the design')
    return(simplify2array(des$model$fisherI[idcs]))
}


#' Update Design
#'
#' \code{update.design} updates the underlying model to ensure the existence of Fisher information matrices for each point in the design.
#' It also updates the corresponding sensitivities.
#'
#' @param des some design.
#'
#' @return \code{update.design} returns an object of \code{class} \code{"design"}. See \code{FedorovWynn} for its structural definition.
#'
#' @seealso \code{\link{reduce}}, \code{\link{update.parm}}, \code{\link{FedorovWynn}}
#'
#' @examples #TODO
#'
#' @export
update.design = function(des) {
    r = des

    mod = update(des$model, des$x)
    r$model = mod

    m = getm(r)
    n = dim(m)[1]
    Mi = solve(getM(m, des$w))

    idcs = getIdcs(des$args$FedorovWynn$names, mod)
    if (identical(idcs, seqi(1, n)))
        sens = sensD(m, Mi)
    else {
        A = getA(idcs, n)
        sens = sensDs(m, Mi, A)
    }

    r$sens = sens
    return(r)
}


# from http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
    if (missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))
}


#' Plot Design
#'
#' \code{plot_design} visualizes the weights and sensitivities of some design.
#' Designs with more than one dimension are projected to a specified margin.
#'
#' The diameter of each circle is proportional to its weight.
#'
#' @param des some design.
#' @param ... other arguments passed to plot.
#' @param margins a vector of indices.
#' Defaults to \code{1}.
#' @param wDes a design from which to take the weights.
#' Defaults to argument \code{des}.
#' See \code{reduce}.
#' @param plus add plus symbols to the sensitivity curve.
#' @param circles draw weights as circles instead of as bars.
#' @param border \code{c(bottom, left, top, right)}, relative margins to add if drawing circles.
#' @param sensArgs a list of arguments for drawing the sensitivity.
#' @param wArgs a list of arguments for drawing the weights.
#'
#' @seealso \code{\link{FedorovWynn}}, \code{\link{reduce}}
#'
#' @examples #TODO
#'
#' @export
plot_design = function(des, ..., margins=NULL, wDes=NULL, plus=T, circles=F, border=c(0.1, 0.1, 0, 0.1), sensArgs=list(), wArgs=list()) {
    if (is.null(margins))
        margins = 1
        #margins = 1:(des$model$dDim)
    if (is.null(wDes))
        wDes = des

    if (1 < length(margins))
        stop('not yet implemented')

    # marginal projections
    x = des$x
    sens = des$sens
    idcs = split(seqi(1, nrow(x)), lapply(margins, function(margin) x[, margin]), drop=T)
    x = x[sapply(idcs, function(idcs) idcs[1]), margins, drop=F]
    sens = sapply(idcs, function(idcs) max(sens[idcs]))

    ord = orderMatrix(x)
    x = x[ord,, drop=F]
    sens = sens[ord]

    wx = wDes$x
    ww = wDes$w
    wsens = wDes$sens
    idcs = split(seqi(1, nrow(wx)), lapply(margins, function(margin) wx[, margin]), drop=T)
    wx = wx[sapply(idcs, function(idcs) idcs[1]), margins, drop=F]
    ww = sapply(idcs, function(idcs) sum(ww[idcs]))
    wsens = sapply(idcs, function(idcs) max(wsens[idcs]))

    args = list(...)
    if (length(margins) == 1) {
        xlim = range(x)
        if (isTRUE(circles)) {
            d = diff(xlim)
            xlim = xlim + c(-1, 1)*d*border[c(2, 4)]
        }
        if (is.null(args$ylim)) {
            ymax = max(sens)
            ylim = c(0, ymax)
            if (isTRUE(circles)) {
                d = diff(ylim)
                ylim = ylim + c(-1, 1)*d*border[c(1, 3)]
            }
        } else {
            ymax = args$ylim[2]
            ylim = NULL
        }
        xlab = colnames(x)
        if (is.null(xlab))
            xlab = paste('x[, c(', toString(margins), ')]', sep='')
        ylab = 'sensitivity'
        args = modifyList(list(xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab), args)

        sensArgs = modifyList(list(col='black'), sensArgs)
        wArgs = modifyList(list(col='black'), wArgs)

        par(mar=c(5, 4, 4, 4) + 0.1)
        do.call(plot, modifyList(list(NA), args))
        do.call(lines, modifyList(list(x, sens), sensArgs))

        if (isTRUE(circles)) {
            do.call(symbols, modifyList(list(wx, rep(0, nrow(wx)), ww, inches=1/2, add=T), wArgs))
            alpha = ww/max(ww)
            idcs = which(1/256 < alpha)
            do.call(abline, modifyList(modifyList(list(v=wx[idcs], lty=3), sensArgs), list(col=add.alpha(sensArgs$col, alpha[idcs])) ))
            #points(wx, rep(0, nrow(wx)), cex=1/2)
        } else {
            if (plus) {
                alpha = ww/max(ww)
                idcs = which(1/256 < alpha)
                wx_ = wx[idcs]
                wsens_ = wsens[idcs]
                if (anyNA(wsens_)) {
                    if (nrow(x) < 2)
                        warning('need at least two design points to interpolate sensitivity')
                    else {
                        nas = is.na(wsens_)
                        wsens_[nas] = approx(x, sens, wx_[nas])$y
                    }
                }
                do.call(points, modifyList(modifyList(list(wx_, wsens_, pch='+'), sensArgs), list(col=add.alpha(sensArgs$col, alpha[idcs])) ))
            }
            par(new=T)
            do.call(plot, modifyList(list(wx, ww, type='h', xlim=args$xlim, ylim=c(0, 1), axes=F, xlab='', ylab=''), wArgs))
            do.call(axis, modifyList(list(4, col.axis=wArgs$col), wArgs))
            do.call(mtext, modifyList(list('weight', side=4, line=3), wArgs))
        }
    }
}


#' D Efficiency
#'
#' \code{Defficiency} computes the D or Ds efficiency measure for some design with respect to some reference design.
#'
#' D efficiency is defined as
#' \deqn{\left(\frac{\left|M(\xi,\bar{\theta})\right|}{\left|M(\xi^{*},\bar{\theta})\right|}\right)^{1/n}}{( abs(M(\xi, \theta))  /  abs(M(\xi*, \theta)) )**(1/n)}
#' and Ds efficiency as
#' \deqn{\left(\frac{\left|M_{11}(\xi,\bar{\theta})-M_{12}(\xi,\bar{\theta})M_{22}^{-1}(\xi,\bar{\theta})M_{12}^{T}(\xi,\bar{\theta})\right|}{\left|M_{11}(\xi^{*},\bar{\theta})-M_{12}(\xi^{*},\bar{\theta})M_{22}^{-1}(\xi^{*},\bar{\theta})M_{12}^{T}(\xi^{*},\bar{\theta})\right|}\right)^{1/s}}{( abs(M11(\xi, \theta) - M12(\xi, \theta) \%*\% solve(M22(\xi, \theta)) \%*\% t(M12(\xi, \theta)))  /  abs(M11(\xi*, \theta) - M12(\xi*, \theta) \%*\% solve(M22(\xi*, \theta)) \%*\% t(M12(\xi*, \theta))) )**(1/s)}
#'
#' @param des some design.
#' @param ref some other design to use as reference.
#'
#' @return \code{Defficiency} returns a single numeric.
#'
#' @seealso \code{\link{FedorovWynn}}, \code{\link{update.design}}
#'
#' @examples #TODO
#'
#' @export
Defficiency = function(des, ref) {
    # TODO check if designs are compatible
    m = getm(des)
    M = getM(m, des$w)

    m = getm(ref)
    n = dim(m)[1]
    Mref = getM(m, ref$w)

    idcs = getIdcs(ref$args$FedorovWynn$names, ref$model)
    if ( identical(idcs, seqi(1, n)) ) {
        print('D')
        return((det(M) / det(Mref)) ** (1/n))
    }

    s = length(idcs)
    M12 = M[idcs, -idcs]
    Mref12 = Mref[idcs, -idcs]

    print('Ds')
    return(( det(M[idcs, idcs] - M12 %*% solve(M[-idcs, -idcs]) %*% t(M12)) / det(Mref[idcs, idcs] - Mref12 %*% solve(Mref[-idcs, -idcs]) %*% t(Mref12)) ) ** (1/s))
}

