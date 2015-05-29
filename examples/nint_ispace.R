library(docopulae)

s = nint_space(nint_scatDim(c(1, 3, 5)),
               nint_gridDim(1:3),
               nint_scatDim(c(2, 4, 6)),
               nint_gridDim(4:5),
               nint_funcDim(function(x) x[1:4]),
               nint_intvDim(-Inf, 3),
               nint_gridDim(6:9),
               nint_intvDim(-Inf, Inf),
               nint_intvDim(-pi, Inf))
s

nint_ispace(s)
# i: indices
# g: data
