

#' Sequence Generation
#'
#' Similar to \code{base::seq}, however \code{by} is strictly 1 by default (not changing with \code{from}, \code{to}) and every closed interval [result, result + \code{by}] is contained in [\code{from}, \code{to}].
#'
#' @param from,to,by see \code{\link[base]{seq}} in \pkg{base}.
#'
#' @return \code{seqi} returns either \code{integer(0)} if no closed interval is contained in [\code{from}, \code{to}] (may be empty) or what an appropriate call to \code{base::seq} returned otherwise.
#'
#' See examples below.
#'
#' @example examples/seq.R
#'
#' @seealso \code{\link[base]{seq}} in \pkg{base}
#'
#' @export
seqi = function(from, to, by=1) {
    length.out_ = as.integer((to - from + sign(by)) / by)
    if (length.out_ <= 0)
        return(integer(0))
    by_ = by
    return(base::seq(from, by=by_, length.out=length.out_))
}


Derivf = function(f, names) {
    temp = Deriv::drule[['[[']]
    assign('[[', list(0), envir=Deriv::drule)

    r = lapply(names, function(name) Deriv::Deriv(f, name))
    names(r) = names

    assign('[[', temp, envir=Deriv::drule)
    return(r)
}

Deriv2f = function(f, names) {
    temp = Deriv::drule[['[[']]
    assign('[[', list(0), envir=Deriv::drule)

    r = replicate(length(names), list())
    base::names(r) = names

    for (i in seq(names)) {
        a = names[[i]]
        print(a)
        d = Deriv::Deriv(f, a)
        for (j in i:length(names)) {
            b = names[[j]]
            print(c(a, b))
            d2 = Deriv::Deriv(d, b)
            r[[a]][[b]] = d2
            r[[b]][[a]] = d2
        }
    }

    assign('[[', temp, envir=Deriv::drule)
    return(r)
}

mirrorMatrix = function(x) {
    r = x
    diag(r) = 0
    r = r + t(r)
    diag(r) = diag(x)
    return(r)
}

is_flat = function(x) !any(sapply(x, inherits, 'list'))

flatten = function(x) {
    if (!inherits(x, 'list'))
        return(list(x))
    if (is_flat(x))
        return(x)
    return(do.call(c, lapply(x, flatten) ))
}

zmin = function(x) ifelse(length(x) == 0, 0, min(x))
zmax = function(x) ifelse(length(x) == 0, 0, max(x))

lproduct = function(x) {
    if (length(x) == 0)
        return(list())
    idcs = lapply(x, function(x) 1:length(x))
    idcs = do.call(expand.grid, idcs)
    colnames(idcs) = names(x[[1]])
    r = apply(idcs, 1, function(idcs) mapply(function(idx, i) x[[i]][[idx]], idcs, 1:length(x), SIMPLIFY=F))
    return(r)
}


#' Integrate Alternative
#'
#' A tolerance wrapper for \code{stats::integrate}. It allows \code{integrate} to reach the maximum number of subdivisions.
#'
#' See \code{stats::integrate}.
#'
#' @seealso \code{\link[stats]{integrate}} in \pkg{stats}
integrateA = function(f, lower, upper, ..., subdivisions=100L, rel.tol=.Machine$double.eps^0.25, abs.tol=rel.tol, stop.on.error=TRUE, keep.xy=FALSE, aux=NULL) {
    r = stats::integrate(f, lower, upper, ..., subdivisions=subdivisions, rel.tol=rel.tol, abs.tol=abs.tol, stop.on.error=F, keep.xy=keep.xy, aux=aux)
    if ( !(r$message %in% c('OK', 'maximum number of subdivisions reached')) ) {
        if (stop.on.error) {
            stop(r$message)
        }
    }
    return(r)
}
