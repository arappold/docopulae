#' Design of Experiments with Copulas
#'
#' A direct approach to optimal designs for copula models based on the Fisher information.
#' Provides flexible functions for building joint PDFs, evaluating the Fisher information and finding optimal designs.
#' It includes an extensible solution to summation and integration called `nint`, functions for transforming, plotting and comparing designs, as well as a set of tools for common low-level tasks.
#'
#' This package builds upon the theoretical result on optimal designs for copula models developed by Elisa Perrone and Werner G. Müller.
#' In their paper named 'Optimal designs for copula models' they introduce an equivalence theorem of Kiefer-Wolfowitz type for D-optimality along with examples and the proof.
#' The proof for D_A-optimality is analogous and is mentioned in an upcoming paper currently under double blind review.
#'
#' @references E. Perrone & W.G. Müller (2016) Optimal designs for copula models, Statistics, 50:4, 917-929, DOI: 10.1080/02331888.2015.1111892
#'
#' @seealso \code{\link{Dsensitivity}}
#'
#' @import graphics grDevices stats utils
#' @importFrom methods substituteDirect
#'
#' @docType package
#' @encoding UTF-8
#' @name docopulae
NULL
