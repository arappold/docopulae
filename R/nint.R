## ideas
# - implement simplify/merge


#' Dimension Type Attribute Values
#'
#' A dimension object is identified by its dimension type attribute \code{"nint_dtype"}.
#' On creation it is set to one of the following.
#' See dimension types in "See Also" below.
#'
#' @format integer
#' @seealso \code{\link{nint_scatDim}}, \code{\link{nint_gridDim}}, \code{\link{nint_intvDim}}, \code{\link{nint_funcDim}}, \code{\link{nint_space}}
#'
#' @name nint_TYPE
NULL

#' @rdname nint_TYPE
#'
#' @usage nint_TYPE_SCAT_DIM = 1
#'
#' @export
nint_TYPE_SCAT_DIM = 1

#' @rdname nint_TYPE
#'
#' @usage nint_TYPE_GRID_DIM = 2
#'
#' @export
nint_TYPE_GRID_DIM = 2

#' @rdname nint_TYPE
#'
#' @usage nint_TYPE_INTV_DIM = 3
#'
#' @export
nint_TYPE_INTV_DIM = 3

#' @rdname nint_TYPE
#'
#' @usage nint_TYPE_FUNC_DIM = 4
#'
#' @export
nint_TYPE_FUNC_DIM = 4


#' Scattered Dimension
#'
#' \code{nint_scatDim} is defined by a vector of values which in combination with other scatter dimensions creates a sparse grid.
#'
#' Imagine using cbind to create a matrix where each row defines a point.
#' This is the role \code{nint_scatDim} plays.
#'
#' @param x a vector of any type.
#'
#' @return \code{nint_scatDim} returns its argument with the dimension type attribute set to \code{nint_TYPE_SCAT_DIM}.
#'
#' @seealso \code{\link{nint_TYPE}}, \code{\link{nint_space}}
#'
#' @export
nint_scatDim = function(x) {
    x = as.vector(x)
    attr(x, 'nint_dtype') = nint_TYPE_SCAT_DIM
    return(x)
}

#' Grid Dimension
#'
#' \code{nint_gridDim} is defined by a vector of values which in combination with other grid dimensions creates a dense grid.
#'
#' Imagine using expand.grid to create a matrix where each row defines a point.
#' This is the role \code{nint_gridDim} plays.
#'
#' @param x a vector of any type.
#'
#' @return \code{nint_scatDim} returns its argument with the dimension type attribute set to \code{nint_TYPE_GRID_DIM}.
#'
#' @seealso \code{\link{nint_TYPE}}, \code{\link{nint_space}}
#'
#' @export
nint_gridDim = function(x) {
    x = as.vector(x)
    attr(x, 'nint_dtype') = nint_TYPE_GRID_DIM
    return(x)
}

#' Interval Dimension
#'
#' \code{nint_intvDim} defines a continous range with a fixed lower and upper limit.
#' The limits may be (negative) \code{Inf}.
#'
#' @param x either a vector containing lower and upper limit, or a single value specifying only the lower limit.
#' @param b the upper limit if \code{x} is the lower limit.
#'
#' @return \code{nint_intvDim} returns a vector of length 2 with the dimension type attribute set to \code{nint_TYPE_INTV_DIM}.
#'
#' @seealso \code{\link{nint_TYPE}}, \code{\link{nint_space}}, \code{\link{nint_transform}}
#'
#' @export
nint_intvDim = function(x, b=NULL) {
    if (is.null(b)) {
        x = as.vector(x)
    } else {
        x = c(x, b)
    }
    if (length(x) != 2) stop('exactly two values shall be specified as limits')
    attr(x, 'nint_dtype') = nint_TYPE_INTV_DIM
    return(x)
}

#' Function Dimension
#'
#' \code{nint_funcDim} defines a functionally dependent dimension.
#' It shall depend solely on the previous dimensions.
#'
#' @param x \code{function(y)}, where \code{y} is the partially realized point in the space. It shall return an object of type \code{nint_intvDim} or a vector.
#'
#' @details Obviously if \code{x} returns an object of type \code{nint_intvDim} the dimension is continous, and discrete otherwise.
#'
#' As the argument of \code{x} contains values for all dimensions the user has to make sure that the function solely depends on values up to the right dimension.
#'
#' @return \code{nint_scatDim} returns its argument with the dimension type attribute set to \code{nint_TYPE_FUNC_DIM}.
#'
#' @seealso \code{\link{nint_TYPE}}, \code{\link{nint_space}}
#'
#' @export
nint_funcDim = function(x) {
    x = as.function(x)
    attr(x, 'nint_dtype') = nint_TYPE_FUNC_DIM
    return(x)
}


nint_dtype = function(x) attr(x, 'nint_dtype')

nint_toString_scatDim = function(x) paste('s(', toString(x), ')', sep='')
nint_toString_gridDim = toString
nint_toString_intvDim = function(x) paste('[', x[1], ', ', x[2], ']', sep='')
nint_toString_funcDim = function(x) paste(format.default(x), sep='', collapse='')

nint_toString_dim = function(x) {
    type = nint_dtype(x)
    if (is.null(type))
        return('<unknown dim type>')
    return(switch(type, nint_toString_scatDim, nint_toString_gridDim, nint_toString_intvDim, nint_toString_funcDim)(x))
}

nint_print_spaceDims = function(x) {
    x = flatten(x)
    if (length(x) == 0)
        return()
    cat(nint_toString_dim(x[[1]]))
    for (i in seqi(2, length(x))) {
        cat(' or ')
        cat(nint_toString_dim(x[[i]]))
    }
    cat('\n')
}


#' Space Validation Errors
#'
#' Error codes for space validation.
#'
#' @format integer
#' @seealso \code{\link{nint_validateSpace}}
#'
#' @name nint_ERROR
NULL

#' @rdname nint_ERROR
#'
#' @usage nint_ERROR_DIM_TYPE = -1001
#'
#' @details \code{nint_ERROR_DIM_TYPE}: dimension type attribute does not exist or is not valid.
nint_ERROR_DIM_TYPE = -1001

#' @rdname nint_ERROR
#'
#' @usage nint_ERROR_SCATTER_LENGTH = -1002
#'
#' @details \code{nint_ERROR_SCATTER_LENGTH}: scatter dimensions of different lengths.
nint_ERROR_SCATTER_LENGTH = -1002

#' @rdname nint_ERROR
#'
#' @usage nint_ERROR_SPACE_TYPE = -1003
#'
#' @details \code{nint_ERROR_SPACE_TYPE}: object not of type \code{"nint_space"}.
nint_ERROR_SPACE_TYPE = -1003

#' @rdname nint_ERROR
#'
#' @usage nint_ERROR_SPACE_DIM = -1004
#'
#' @details \code{nint_ERROR_SPACE_DIM}: subspaces with different number of dimensions.
nint_ERROR_SPACE_DIM = -1004


nint_validateSpaceDim = function(x, refl) {
    # check - dim type
    #       - scat length
    type = nint_dtype(x)
    if (is.null(type))
        return(nint_ERROR_DIM_TYPE)
    pass = function(x, refl) 0
    f = switch(type, function(x, refl) ifelse(length(x) != refl, nint_ERROR_SCATTER_LENGTH, 0), pass, pass, pass)
    if (is.null(f))
        return(nint_ERROR_DIM_TYPE)
    return(f(x, refl))
}

nint_validateSpaceDims = function(x, refl) zmin(sapply(flatten(x), nint_validateSpaceDim, refl))


#' Space
#'
#' \code{nint_space} defines an n-dimensional space as a list of dimensions.
#' A space may contain subspaces.
#' A space without subspaces is called true subspace.
#'
#' If a space contains at least one list structure of dimension objects it consists of subspaces each defined by a combination of dimension objects along all dimensions.
#' See \code{\link{nint_expandSpace}} on how to extract them.
#'
#' @param ... dimensions each of which may be an actual dimension object or a list structure of dimension objects.
#'
#' @return \code{nint_space} returns an object of \code{class} \code{"nint_space"}.
#'
#' @seealso \code{\link{nint_scatDim}}, \code{\link{nint_gridDim}}, \code{\link{nint_intvDim}}, \code{\link{nint_funcDim}}, \code{\link{nint_transform}}, \code{\link{nint_integrate}}, \code{\link{nint_expandSpace}}
#'
#' @example examples/nint_space.R
#'
#' @export
nint_space = function(...) {
    r = list(...)
    class(r) = 'nint_space'
    return(r)
}

#' Print Space
#'
#' Prints object of type \code{nint_space} in a convenient way.
#'
#' Each line represents a dimension.
#' Format: "<dim idx>: <dim repr>".
#' Each dimension has its own representation which should be easy to understand.
#' \code{nint_scatDim} representations are marked by \code{"s()"}.
#'
#' @param x some space.
#' @param ... ignored.
#'
#' @seealso \code{\link{nint_space}}
#'
#' @export
print.nint_space = function(x, ...) {
    for (i in seqi(1, length(x))) {
        cat(i, ': ', sep='')
        nint_print_spaceDims(x[[i]])
    }
}

#' Validate Space
#'
#' \code{nint_validateSpace} performs some checks to make sure the space is well defined.
#'
#' @param x some space.
#'
#' @return \code{nint_validateSpace} returns 0 if everything is fine, or an error code.
#' See \code{\link{nint_ERROR}}.
#'
#' @seealso \code{\link{nint_ERROR}}, \code{\link{nint_space}}
#'
#' @example examples/nint_validateSpace.R
#'
#' @export
nint_validateSpace = function(x) {
    x = flatten(x)
    if (length(x) == 0)
        return(0)
    return(min(sapply(x, nint_validateSpace_, x[[1]])))
}

nint_validateSpace_ = function(x, refs=NULL) {
    if (!inherits(x, 'nint_space')) # check class
        return(nint_ERROR_SPACE_TYPE)
    if (length(x) != length(refs)) # check nint_space dimension
        return(nint_ERROR_SPACE_DIM)
    x = lapply(x, flatten)
    return(zmin( sapply(x, nint_validateSpaceDims, nint_validateSpace_getRefl(x)) ))
}

nint_validateSpace_getRefl = function(x) {
    return(zmax(sapply(x, function(x) zmax(sapply(x, nint_validateSpace_getRefl_) )) ))
}

nint_validateSpace_getRefl_ = function(x) {
    type = nint_dtype(x)
    if (is.null(type) || type != nint_TYPE_SCAT_DIM)
        return(0)
    return(length(x))
}

#' Expand Space
#'
#' \code{nint_expandSpace} expands some space or list structure of spaces to a list of all true subspaces.
#'
#' @param x some space or a list structure of spaces.
#'
#' @return \code{nint_expandSpace} returns a list of objects each of \code{class} \code{"nint_space"}.
#' Each space is a true subspace.
#'
#' @seealso \code{\link{nint_space}}
#'
#' @example examples/nint_expandSpace.R
#'
#' @export
nint_expandSpace = function(x) {
    return(flatten(lapply(flatten(x), nint_expandSpace_) ))
}

nint_expandSpace_ = function(x) {
    return( lapply(lproduct(lapply(x, flatten)), function(x) do.call(nint_space, x)) )
}


## Create Integration Space
##
## \code{nint_ispace} transforms a true subspace into a data structure which is used by \code{nint_integrate} to efficiently integrate over the entire space.
##
## Assumes interchangeability of dimensions (except function dimensions).
##
## @param x a true subspace.
##
## @return \code{nint_ispace} returns a list.
##
## @seealso \code{\link{nint_space}}, \code{\link{nint_integrate}}
##
## @example examples/nint_ispace.R
##
## @export
nint_ispace = function(x) {
    # - group dimensions by type and put them in the following order
    #   - scatter (s)
    #   - grid (g)
    #   - interval (i)
    #   - function (f)
    # - bind scattered dimensions by column
    # - expand grid dimensions to grid
    # - bind interval limits by row
    # - create function list
    #
    # - result := list(type=list(i=idcs, g=data))
    #   type := c('s', 'g', 'i', 'f')
    #   data := matrix      if type == 's'
    #        or data.frame  if type == 'g' or type == 'i'
    #        or list        if type == 'f'

    if (length(x) == 0)
        return(list())

    rr = list()
    si = integer(0) # scatter
    gi = integer(0) # grid
    ii = integer(0) # interval
    fi = integer(0) # function

    for (i in seqi(1, length(x))) {
        xx = x[[i]]
        if (inherits(xx, 'list'))
            stop('argument is no true subspace')
        type = nint_dtype(xx)
        if (is.null(type))
            stop('unknown dimension type')

        if (type == nint_TYPE_FUNC_DIM) {
            if (length(si) != 0) rr = c(rr, list(s=si))
            if (length(gi) != 0) rr = c(rr, list(g=gi))
            if (length(ii) != 0) rr = c(rr, list(i=ii))
            si = integer(0)
            gi = integer(0)
            ii = integer(0)
            fi = c(fi, i)
            next
        }
        if (length(fi) != 0) {
            rr = c(rr, list(f=fi))
            fi = integer(0)
        }

        if (type == nint_TYPE_SCAT_DIM) {
            si = c(si, i)
        } else if (type == nint_TYPE_GRID_DIM) {
            gi = c(gi, i)
        } else if (type == nint_TYPE_INTV_DIM) {
            ii = c(ii, i)
        } else {
            stop('unknown dim type')
        }
    }

    if (length(si) != 0) rr = c(rr, list(s=si))
    if (length(gi) != 0) rr = c(rr, list(g=gi))
    if (length(ii) != 0) rr = c(rr, list(i=ii))
    if (length(fi) != 0) rr = c(rr, list(f=fi))

    ts = function(i) do.call(cbind, x[i])
    tg = function(i) do.call(expand.grid, x[i])
    ti = function(i) do.call(rbind, x[i])
    tf = function(i) x[i]
    tt = list(s=ts, g=tg, i=ti, f=tf)
    r = mapply(function(type, i) list(i=i, g=tt[[type]](i)), names(rr), rr, SIMPLIFY=F)

    ## create matrix containing the depth (first column) and position (second column) for each dimension
    #depth = rep.int(1:length(rr), sapply(rr, length))[order(unlist(rr))]
    #attr(r, 'depth') = cbind(depth, mapply(match, 1:length(x), lapply(depth, function(d) r[[d]]$i)), deparse.level=0)
    return(r)
}

nint_ispaces = function(x) {
    r = nint_validateSpace(x)
    if (r != 0)
        stop(r)
    x = nint_expandSpace(x)
    return(lapply(x, nint_ispace))
}


# replaced by nint_transform
### Transform Infinite Intervals
###
### \code{nint_transformInf} transforms (semi) infinite intervals to finite intervals.
###
### Transformations: \itemize{
### \item (\code{-Inf}, a): \code{x = a + 1 - 1/t} with t between 0 and 1
### \item (a, \code{Inf}): \code{x = a - 1 + 1/t} with t between 0 and 1
### \item (\code{-Inf}, \code{Inf}): \code{x = t/(1 - abs(t))} with t between -1 and 1
### }
###
### @param f a scalar-valued function.
### @param ispace an integration space.
###
### @return \code{nint_transformInf} returns a list containing the transformed function and transformed integration space.
###
### @seealso \code{\link{nint_integrate}}, \code{\link{nint_ispace}}
###
### @example examples/nint_transformInf.R
###
### @export
#nint_transformInf = function(f, ispace) {
    #r = list()
    #rispace = ispace # copy

    #types = names(ispace)
    #for (i in seq(ispace)) {
        #type = types[[i]]
        #if (type != 'i') {
            #next
        #}

        #stage = ispace[[i]]

        ## __  .. finite limits
        ## l  .. -Inf, a
        ## u  .. a, Inf
        ## b .. -Inf, Inf
        #infTypes = apply(stage$g, 1, function(limits) ifelse(limits[[1]] == -Inf, ifelse(limits[[2]] == Inf, 'b', 'l'), ifelse(limits[[2]] == Inf, 'u', '__')))

        #lu = rep(F, length(infTypes))
        #j = infTypes == 'l'
        #lu = lu | j
        #if (any(j))
            #r$l = rbind(r$l, cbind(stage$i[j], stage$g[j, 2]))
        #j = infTypes == 'u'
        #lu = lu | j
        #if (any(j))
            #r$u = rbind(r$u, cbind(stage$i[j], stage$g[j, 1]))
        #if (any(lu))
            #rispace[[i]]$g[lu,] = matrix(c(0, 1), nrow=sum(lu), ncol=2, byrow=T)

        #j = infTypes == 'b'
        #if (any(j)) {
            #r$b = c(r$b, stage$i[j])
            #rispace[[i]]$g[j,] = matrix(c(-1, 1), nrow=sum(j), ncol=2, byrow=T)
        #}
    #}

    #if (!is.null(r$l) && nrow(r$l) != 0) {
        #li = r$l[,1] # idcs
        #ll = r$l[,2] # upper limits
        #lf = f
        #f = function(x, ...) {
            #t1 = 1 / x[li]
            #t2 = prod(t1 * t1)
            #if (t2 == Inf)
                #return(0)
            #x[li] = ll + 1 - t1
            #return(lf(x, ...) * t2)
        #}
    #}
    #if (!is.null(r$u) && nrow(r$u) != 0) {
        #ui = r$u[,1] # idcs
        #ul = r$u[,2] # lower limits
        #uf = f
        #f = function(x, ...) {
            #t1 = 1 / x[ui]
            #t2 = prod(t1 * t1)
            #if (t2 == Inf)
                #return(0)
            #x[ui] = ul - 1 + t1
            #return(uf(x, ...) * t2)
        #}
    #}
    #if (!is.null(r$b)) {
        #bi = r$b
        #bf = f
        #f = function(x, ...) {
            #t1 = x[bi]
            #t2 = 1 / (1 - abs(t1))
            #t3 = prod(t2 * t2)
            #if (t3 == Inf)
                #return(0)
            #x[bi] = t1 * t2
            #return(bf(x, ...) * t3)
        #}
    #}

    #return(list(f=f, ispace=rispace))
#}


# not needed
## @export
#nint_applyToSpace = function(space, d, f) {
    #r = lapply(flatten(space), nint_applyToSpace_, d, f)
    #r = flatten(r)
    #if (length(r) == 1)
        #return(r[[1]])
    #return(r)
#}

#insert = function(x, idcs, l) {
    ## x .. list of objects
    ## idcs .. vector of indices
    ## l .. list to insert to
    #r = l
    #d = length(x) - length(idcs)
    #if (d == 0) {
        #r[idcs] = x
    #} else if (d < 0) { # more idcs than items
        #r[idcs[seqi(1, length(x))]] = x # set
        #r = r[-idcs[seqi(length(x) + 1, length(idcs))]] # del
        #r = do.call(nint_space, r)
    #} else { # more items than idcs
        #r[idcs] = x[seqi(1, length(idcs))] # set
        #last = idcs[length(idcs)]
        #r = c(r[seqi(1, last)], x[seqi(length(idcs) + 1, length(x))], r[seqi(last + 1, length(r))]) # insert
        #r = do.call(nint_space, r)
    #}
    #return(r)
#}

#nint_applyToSpace_ = function(space, d, f) {
    #dds = lapply(space[d], flatten) # list(list(dim))
    #r = lapply(lproduct(dds), f) # list(obj)
    #r = lapply(r, function(x) if (inherits(x, 'list')) x else list(x)) # list(list(dim))
    #r = lapply(r, insert, d, space) # list(space)
    #return(r)
#}


ratiog = function(x) {
    s = sign(x)
    if (any(s == 0, na.rm=T))
        s[s == 0] = 1
    r = ifelse(is.infinite(x), s, x/abs(x + s))
}

transforms = list(tan=list(g=atan, gi=tan, gij=function(x) 1 + tan(x)**2),
                  ratio=list(g=ratiog, gi=function(x) x / (1 - abs(x)), gij=function(x) 1 / (1 - abs(x))**2))

#' Transform Integral
#'
#' \code{nint_transform} applies monotonic continous transformations.
#' A common purpose is to transform infinite intervals to finite ones.
#'
#' If the transformation is vector valued, that is \code{y = c(y1, ..., yn) = g(c(x1, ..., xn))}, then each component of \code{y} shall exclusively depend on the corresponding component of \code{x}.
#' An incorrect expression for this would be: \code{y[i] = g[i](x[i])}.
#'
#' Builtins: \itemize{
#' \item tan: \code{y = atan(x)}
#' \item ratio: \code{y = x/abs(x + sign(x))} with \code{sign(0) == 1}
#' }
#'
#' @param f \code{function(x, ...)}.
#' @param space some space.
#' @param dIdcs an integer vector of indices, the dimensions to transform.
#' @param trans either a name of some builtin transformation or \code{list(g=function(x), gi=function(y), gij=function(y))} where \code{y = g(x)}, \code{gi(y) = x} and \code{gij} is the first derivative of \code{gi} with respect to \code{y}.
#' @param infZero the value to return if the jacobian is infinite and \code{f} returns \code{0}.
#'
#' @return \code{nint_transform} returns a named list containing the transformed function and space.
#'
#' @seealso \code{\link{nint_space}}, \code{\link{nint_integrate}}, \code{\link{nint_intvDim}}
#'
#' @example examples/nint_transform.R
#'
#' @export
nint_transform = function(f, space, dIdcs, trans, infZero=0) {
    tt = list(f, space, dIdcs, trans, infZero)

    if (is.character(trans)) {
        tt = transforms[[trans]]
        if (is.null(tt))
            stop(paste('unknown transformation \'', trans, '\'', sep=''))
        trans = tt
    }

    g = trans[['g']]
    gi = trans[['gi']]
    gij = trans[['gij']]
    neg = 0

    rspace = lapply(flatten(space), function(space) { # for each true space
        x = rep(NA, length(dIdcs))
        rr = mapply(function(i, d) { # for each dimension
            rapply(list(d), function(d) { # for each true dimension
                dtype = attr(d, 'nint_dtype')
                if (is.null(dtype) || dtype != nint_TYPE_INTV_DIM)
                    stop('dimensions shall exclusively be of type interval')

                v = sapply(d, function(xx) { # transform values
                    x[i] = xx
                    return(g(x)[i])
                })

                if (diff(v) < 0) {
                    neg <<- neg + 1
                    v = v[c(2, 1)]
                }
                return(nint_intvDim(v))
            }, how='replace')[[1]]
        }, seqi(1, length(dIdcs)), space[dIdcs], SIMPLIFY=F)
        r = space
        r[dIdcs] = rr
        return(r)
    })

    if (length(rspace) == 1)
        rspace = rspace[[1]]

    sig = ifelse(neg %% 2, -1, 1)
    rf = function(x, ...) {
        xx = x[dIdcs]
        x[dIdcs] = gi(xx)
        v = f(x, ...)
        j = prod(gij(xx))
        if (is.infinite(j)) {
            if (v == 0) {
                v = infZero
                j = sign(j)
            }
        }
        return(sig*v*j)
    }

    return(list(f=rf, space=rspace))
}



#' Integrate N Cube
#'
#' Interface to the integration over interval dimensions.
#'
#' @usage nint_integrateNCube(f, lowerLimit, upperLimit, ...)
#'
#' @details \code{nint_integrateNCube} is a reference to the function that \code{nint_integrate} calls with exactly these arguments to integrate over interval dimensions.
#' See examples below on how to replace it with a different function.
#'
#' @param f the scalar-valued function (integrand) to be integrated.
#' @param lowerLimit the lower limits of integration.
#' @param upperLimit the upper limits of integration.
#' @param ... other arguments passed to \code{f}.
#'
#' @return \code{nint_integrateNCube} returns a single numeric.
#'
#' @seealso \code{\link{nint_integrate}}
#'
#' @example examples/nint_integrateNCube.R
#'
#' @name nint_integrateNCube
#'
#' @export
NULL

#' @param integrate \code{function(f, lowerLimit, upperLimit, ...)} which calls \code{stats::integrate}.
#'
#' @details \code{nint_integrateNCube_integrate} uses \code{integrate} recursively.
#' Downside: number of function evaluations is \code{(subdivisions * 21) ** N}.
#' This is the default because no package is required.
#' However, you most likely want to consider different solutions.
#'
#' @seealso \code{\link{integrateA}}, \code{\link[stats]{integrate}} in \pkg{stats}
#'
#' @rdname nint_integrateNCube
#'
#' @export
nint_integrateNCube_integrate = function(integrate) {
    tt = integrate
    r = function(f, lowerLimit, upperLimit, ...) {
        maxDepth = length(lowerLimit)

        y = rep(0, maxDepth)
        d = 0
        g = Vectorize(function(x, ...) {
            y[d] <<- x
            if (d == maxDepth) {
                return(f(y, ...))
            }

            d <<- d + 1
            r = integrate(g, lowerLimit[d], upperLimit[d], ...)$value
            d <<- d - 1
            return(r)
        }, 'x')

        return(g(0, ...))
    }
    return(r)
}

#' @param adaptIntegrate \code{function(f, lowerLimit, upperLimit, ...)} which calls \code{cubature::adaptIntegrate}.
#'
#' @details \code{nint_integrateNCube_cubature} is a trivial wrapper for \code{cubature::adaptIntegrate}.
#'
#' @seealso \code{\link[cubature]{adaptIntegrate}} in package \pkg{cubature}
#'
#' @rdname nint_integrateNCube
#'
#' @export
nint_integrateNCube_cubature = function(adaptIntegrate) {
    tt = adaptIntegrate
    r = function(f, lowerLimit, upperLimit, ...) {
        return(adaptIntegrate(f, lowerLimit, upperLimit, ...)$integral)
    }
    return(r)
}

#' @param createIntegrationGrid \code{function(dimension)} which calls \code{SparseGrid::createIntegrationGrid}.
#'
#' @details \code{nint_integrateNCube_SparseGrid} is an almost trivial wrapper for \code{SparseGrid::createIntegrationGrid}.
#' It scales the grid to the integration region.
#'
#' @seealso \code{\link[SparseGrid]{createIntegrationGrid}} in package \pkg{SparseGrid}
#'
#' @rdname nint_integrateNCube
#'
#' @export
nint_integrateNCube_SparseGrid = function(createIntegrationGrid) {
    tt = createIntegrationGrid
    r = function(f, lowerLimit, upperLimit, ...) {
        n = length(lowerLimit)
        d = upperLimit - lowerLimit
        grid = createIntegrationGrid(n) # generates grid on [0, 1] ** n
        # transform to [lowerLimit, upperLimit]
        grid$nodes = sweep(sweep(grid$nodes, 2, d, '*'), 2, lowerLimit, '+')

        r = apply(grid$nodes, 1, f, ...)
        r = (r %*% grid$weights)[1, 1] * prod(d)
        return(r)
    }
    return(r)
}

nint_integrateNCube = nint_integrateNCube_integrate(integrateA)


#' Integrate N Function
#'
#' Inferface to the integration over function dimensions.
#'
#' @usage nint_integrateNFunc(f, funcs, x0, i0, ...)
#'
#' @details \code{nint_integrateNCube} is a reference to the function that \code{nint_integrate} calls with exactly these arguments to integrate over function dimensions.
#' See examples below on how to replace it with a different function.
#'
#' @param f the scalar-valued function (integrand) to be integrated.
#' @param funcs the list of functions. See \code{\link{nint_funcDim}}.
#' @param x0 the partially realized point.
#' @param i0 the vector of indices where the function dimensions fit into \code{x0}.
#' @param ... other arguments passed to \code{f}.
#'
#' @return \code{nint_integrateNFunc} returns a single numeric.
#'
#' @seealso \code{\link{nint_integrate}}
#'
#' @example examples/nint_integrateNFunc.R
#'
#' @name nint_integrateNFunc
#'
#' @export
NULL

#' @param integrateIntv \code{function(f, lowerLimit, upperLimit, ...)}, which performs one dimensional integration.
#'
#' @details \code{nint_integrateNFunc_recursive} returns a recursive implementation which directly sums over discrete dimensions and uses \code{integrateIntv} otherwise. In conjunction with \code{stats::integrate} this is the default.
#'
#' @seealso \code{\link[stats]{integrate}} in \pkg{stats}
#'
#' @rdname nint_integrateNFunc
nint_integrateNFunc_recursive = function(integrateIntv) {
    r = function(f, funcs, x0, i0, ...) {
        maxDepth = length(funcs)
        if (maxDepth == 0)
            return(0)

        y = rep(0, maxDepth)
        d = 0
        g = Vectorize(function(x, ...) {
            y[d] <<- x
            x0[i0[d]] <<- x
            if (d == maxDepth) {
                return(f(y, ...))
            }

            d <<- d + 1
            n = funcs[[d]](x0)
            type = nint_dtype(n)
            if (type == nint_TYPE_INTV_DIM)
                r = integrateIntv(g, n[1], n[2], ...)
            else
                r = sum(vapply(n, g, 0, ...))
            d <<- d - 1
            return(r)
        }, 'x')

        return(g(0, ...))
    }
    return(r)
}

nint_integrateNFunc = nint_integrateNFunc_recursive(function(...) integrateA(...)$value)



#' Integrate
#'
#' \code{nint_integrate} performs n-dimensional integration of some scalar-valued function over the supplied space.
#'
#' \code{nint_integrate} calls \code{nint_integrateNCube} and \code{nint_integrateNFunc} to integrate over interval and function dimensions.
#' See their help pages on how to replace them by different functions.
#'
#' Exploration of the space is done recursively with the order of dimensions optimized for efficiency.
#' Therefore interchangeability of dimensions (except for function dimensions) is assumed.
#'
#' @param f the scalar-valued function (integrand) to be integrated.
#' @param space some space.
#' @param ... other arguments passed to \code{f}.
#'
#' @return \code{nint_integrate} returns a single numeric.
#'
#' @seealso \code{\link{nint_space}}, \code{\link{nint_transform}}, \code{\link{nint_integrateNCube}}, \code{\link{nint_integrateNFunc}}
#'
#' @example examples/nint_integrate.R
#'
#' @export
nint_integrate = function(f, space, ...) {
    ispaces = nint_ispaces(space)
    if (length(ispaces) == 0)
        return(0)

    # globals
    y = rep(0, sum(sapply(ispaces[[1]], function(x) length(x$i)) ))
    ispace = NULL
    gg = list() # sequence of functions
    maxDepth = 0
    d = 0 # depth
    i = 0 # indices
    g = NULL # data

    gs = function(...) sum(apply(g, 1, gd, ...)) # g scatter
    gi = function(...) nint_integrateNCube(gd, g[,1], g[,2], ...) # g interval
    gf = function(...) nint_integrateNFunc(gd, g, y, i, ...) # g function
    gl = list(s=gs, g=gs, i=gi, f=gf) # g list

    gd = function(x, ...) { # g dispatch
        y[i] <<- x
        if (d == maxDepth)
            return(f(y, ...))

        # save state
        ti = i
        tg = g

        # prepare descend
        d <<- d + 1
        idim = ispace[[d]]
        i <<- idim$i
        g <<- idim$g

        # descend
        r = gg[[d]](...)

        # restore state
        d <<- d - 1
        i <<- ti
        g <<- tg
        return(r)
    }

    af = f # argument
    r = 0
    for (ispace in ispaces) {
        gg = gl[names(ispace)]
        maxDepth = length(ispace)

        r = r + gd(0, ...)
    }

    return(r)
}



