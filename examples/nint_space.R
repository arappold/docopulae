library(docopulae)

s = nint_space(nint_gridDim(1:4),
               list(nint_intvDim(-Inf, -3), nint_intvDim(3, Inf)),
               list(list(list(nint_scatDim(c(1, 2))), nint_scatDim(c(2, 3))), nint_scatDim(c(3, 4))))
s
