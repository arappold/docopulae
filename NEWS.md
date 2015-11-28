# docopulae 0.3.3

* fixed defaults for `method.args` for `numDeriv2Logf`
* `plot.desigh` won't draw second axis anymore if `sensArgs$axes == F`

# docopulae 0.3.2

* `numDerivLogf` and `numDeriv2Logf` optionally take a log likelihood function
* sensitivity function redefined to a positive function, concerns
  * `Dsensitivity`
  * argument `tol` to `FedorovWynn`
  * new argument `sensTol` to `plot.desigh`

# docopulae 0.3.1

* `buildf` is more general

# docopulae 0.3.0

* class `desigh` redefined, constructor `design` added
* major change in workflow
  * sensitivity function is 'liberated'
  * `Dsensitivity` builds a sensitivity function for D- and Ds-optimality
  * `FedorovWynn` takes a sensitivity function
  * related functions and their arguments are adjusted accordingly
* component `dDim` removed from model definition
* `update.param` also takes designs therefore replacing `update.desigh` and `update_reference`
* new defaults for arguments `dsNames` and `names`

# docopulae 0.2.0

* initial CRAN submission
* start package versioning
* all previous commits are considered irrelevant
