% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/nint.R
\docType{data}
\name{nint_ERROR}
\alias{nint_ERROR}
\alias{nint_ERROR_DIM_TYPE}
\alias{nint_ERROR_SCATTER_LENGTH}
\alias{nint_ERROR_SPACE_DIM}
\alias{nint_ERROR_SPACE_TYPE}
\title{Space Validation Errors}
\format{integer}
\usage{
nint_ERROR_DIM_TYPE = -1001

nint_ERROR_SCATTER_LENGTH = -1002

nint_ERROR_SPACE_TYPE = -1003

nint_ERROR_SPACE_DIM = -1004
}
\description{
Error codes for space validation.
}
\details{
\code{nint_ERROR_DIM_TYPE}: dimension type attribute does not exist or is not valid.

\code{nint_ERROR_SCATTER_LENGTH}: scatter dimensions of different lengths.

\code{nint_ERROR_SPACE_TYPE}: object not of type \code{"nint_space"}.

\code{nint_ERROR_SPACE_DIM}: subspaces with different number of dimensions.
}
\seealso{
\code{\link{nint_validateSpace}}
}
\keyword{datasets}
