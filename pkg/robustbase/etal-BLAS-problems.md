# Problems reported by CRAN (Brian Ripley) when checking robustbase with non-R versions of BLAS:

- *NOTA BENE* /usr/bin/R on Fedora has "BLAS & LAPACK = "OpenBLAS" :
   `/usr/lib64/R/lib/libRblas.so` (when open in emacs, you can search for `OpenBLAS`)

## Overview (manual, from the ~/R/Pkgs/robustbase/robustbase_<kind>.out  file):

Update after 0.93-3 was released to CRAN, 2018-..-.. :

| BLAS kind | file                  | error msg                                          | error kind        | Ok | comments |
| :-------- | :-------------------- | :------------------------------------------------  | :---------------  | -- | -------- |
| ATLAS     | tests/NAcoef.R        | predict(rm1, se.fit = TRUE, interval....           |                   |    |          |
| ATLAS     | tests/nlregrob-tst.R  | update(Cfit.40out, data=d.exp.Hlev)-> nls err [^1] |                   |    |          |
| MKL       | tests/nlregrob-tst.R  | update(Cfit.40out, data=d.exp.Hlev)-> nls err [^1] |                   |    |          |
| OpenBLAS  | tests/nlregrob-tst.R  | update(Cfit.40out, ==== same as MKL above ....     |                   |    |          |

### Ok: `x`: fixed for 0.93-xx -- `no`: still in 0.9...;  X: disabled for 0.9.....


[^1]: Error in nls(formula = y ~ Expo(x, a, b), data = d.exp.Hlev, start = c(a = 1,  :
    step factor 0.000488281 reduced below 'minFactor' of 0.000976562
