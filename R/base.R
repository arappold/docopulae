

#' Sequence Generation
#'
#' \code{seq1} is similar to \code{base::seq}, however \code{by} is strictly \code{1} by default and \code{integer(0)} is returned if range is empty.
#'
#' @param from,to,by see \code{\link[base]{seq}}.
#'
#' @return \code{seq1} returns either \code{integer(0)} if range is empty or what an appropriate call to \code{base::seq} returns otherwise.
#'
#' See examples below.
#'
#' @example examples/seq1.R
#'
#' @seealso \code{\link[base]{seq}}
#'
#' @export
seq1 = function(from, to, by=1) {
    if (to == from)
        return(from)
    if (sign(to - from) != sign(by))
        return(integer(0))
    return(seq(from, to, by))
}


Derivf = function(f, names) {
    # applies Deriv
    # ignores [[
    #temp = Deriv::drule[['[[']]
    #assign('[[', list(0), envir=Deriv::drule)

    r = lapply(names, function(name) Deriv::Deriv(f, name))
    base::names(r) = names

    #assign('[[', temp, envir=Deriv::drule)
    return(r)
}

Deriv2f = function(f, names) {
    # applies Deriv twice
    # ignores [[
    #temp = Deriv::drule[['[[']]
    #assign('[[', list(0), envir=Deriv::drule)

    r = replicate(length(names), list())
    base::names(r) = names

    for (i in seq(names)) {
        a = names[[i]]
        d = Deriv::Deriv(f, a)
        for (j in i:length(names)) {
            b = names[[j]]
            d2 = Deriv::Deriv(d, b)
            r[[a]][[b]] = d2
            r[[b]][[a]] = d2
        }
    }

    #assign('[[', temp, envir=Deriv::drule)
    return(r)
}


mirrorMatrix = function(x) {
    # transforms upper/lower diagonal matrix to full matrix
    r = x
    diag(r) = 0
    r = r + t(r)
    diag(r) = diag(x)
    return(r)
}


## nested list helper
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
    # product expands a list of lists
    if (length(x) == 0)
        return(list())
    idcs = lapply(x, seq)
    idcs = do.call(expand.grid, idcs)
    colnames(idcs) = names(x[[1]])
    r = apply(idcs, 1, function(idcs)
        mapply(function(idx, i) x[[i]][[idx]], idcs, 1:length(x), SIMPLIFY=F))
    return(r)
}


#' Integrate Alternative
#'
#' \code{integrateA} is a tolerance wrapper for \code{stats::integrate}.
#' It allows \code{integrate} to reach the maximum number of subdivisions.
#'
#' See \code{\link[stats]{integrate}}.
#'
#' @param f,lower,upper,...,subdivisions,rel.tol,abs.tol,stop.on.error,keep.xy,aux see \code{\link[stats]{integrate}}.
#'
#' @seealso \code{\link[stats]{integrate}}
#'
#' @export
integrateA = function(f, lower, upper, ..., subdivisions=100L, rel.tol=.Machine$double.eps^0.25, abs.tol=rel.tol, stop.on.error=TRUE, keep.xy=FALSE, aux=NULL) {
    r = stats::integrate(f, lower, upper, ..., subdivisions=subdivisions, rel.tol=rel.tol, abs.tol=abs.tol, stop.on.error=F, keep.xy=keep.xy, aux=aux)
    if ( !(r$message %in% c('OK', 'maximum number of subdivisions reached')) ) {
        if (stop.on.error) {
            stop(r$message)
        }
    }
    return(r)
}


clusterPeak = function(x, y, maxDist) {
    # x = row matrix of points
    # y = corresponding vector of values (of length nrow(x))
    #
    # iteratively assigns points to the nearest cluster
    # in reverse order of magnitude of y

    r = rep(0, length(y))

    dists = as.matrix(dist(x))
    idcs = seq1(1, length(y))

    rr = 1

    while (length(idcs) != 0) {
        i = which.max(y[idcs])
        idx = idcs[i] # idx of max y
        ds = dists[,idx] # distances to all other points

        cIdcs = which(ds <= maxDist) # candidates
        cIdcs = cIdcs[order(ds[cIdcs])] # sort by distance
        cSub = match(T, r[cIdcs] != 0) # first already matched
        if (is.na(cSub)) { # none matched
            r[idx] = rr
            rr = rr + 1
            next
        }

        r[idx] = r[cIdcs[cSub]] # match
        idcs = idcs[-i]
    }

    return(r)
}


# row order of a matrix
orderMatrix = function(x) do.call(order, lapply(seq1(1, ncol(x)), function(i) x[,i]))

indexMatrix = function(x, y) {
    # for each row in y, finds the index of first matching row in x
    ordx = orderMatrix(x)
    ordy = orderMatrix(y)
    r = which(duplicated(rbind(y[ordy,, drop=F], x[ordx,, drop=F]))) - nrow(y)
    # here r contains indices of matching rows in sorted x
    return(ordx[r][order(ordy)])
}

