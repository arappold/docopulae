## Test environments

* local Linux Mint install, R 3.1.2
* win-builder (devel and release)

## R CMD check results

This is the first submission.

There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking S3 generic/method consistency ... NOTE
  Found the following apparent S3 methods exported but not registered:
    plot.design

  `"design"` as well as `plot(design)` is very natural in the context of this package.
  `docopulae::plot.design` would not 'override' `graphics::plot.design` if it were registered.

## Downstream dependencies

There are currently no downstream dependencies for this package.
