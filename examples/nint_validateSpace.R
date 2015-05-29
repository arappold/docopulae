library(docopulae)

## valid spaces
s = nint_space()
s
nint_validateSpace(s)

s = nint_space(nint_intvDim(-1, 1))
s
nint_validateSpace(s)

## -1001
s = nint_space(1)
s
nint_validateSpace(s)

## -1002
s = nint_space(list(nint_scatDim(c(1, 2)), nint_scatDim(c(1, 2, 3))))
s
nint_validateSpace(s)

s = nint_space(nint_scatDim(c(1, 2)),
               nint_scatDim(c(1, 2, 3)))
s
nint_validateSpace(s)

## -1003
nint_validateSpace(1)
nint_validateSpace(list(nint_space())) # valid
nint_validateSpace(list(1))

## -1004
s1 = nint_space(nint_intvDim(-Inf, Inf),
                nint_intvDim(-3, 2))
s2 = do.call(nint_space, s1[1])
s1
s2
nint_validateSpace(list(s1, s2))
