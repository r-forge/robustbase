# Problems reported by CRAN (Brian Ripley) when checking robustbase with non-R versions of BLAS:

## Overview (manual, from the ~/R/Pkgs/robustbase/robustbase_<kind>.out  file):

Update after 0.93-0 was released to CRAN, 2018-04-22 : The three '?' actually did become _no_

| BLAS kind | file                  | error msg                                           | error kind        | Ok | comments |
| :-------- | :-------------------- | :------------------------------------------------   | :---------------  | -- | -------- |
| ATLAS     | tests/mc-strict.R     | a5$nonOut are not all TRUE                          | diff. outlier set | no | disabled |
| ATLAS     | tests/nlrob-tst.R     | assert.EQ(coef(nlR1),.. tol=1e-9) -- gave 1.073e-9  |                   | x  |          |
| MKL       | tests/nlregrob-tst.R  | update(Cfit.40out, data=d.exp.Hlev)-> nls err [^1] |                    | no |   TODO   |
| OpenBLAS  | tests/m-s-estimator.R | all.equal(mC[nm], mR[nm], tol=5e-15)) -> [^2]       |                   | x  |          |
| OpenBLAS  | tests/nlregrob-tst.R  | update(Cfit.40out, ==== same as MKL above ....      |                   | no |          |

### Ok: `x`: fixed for 0.93-0 -- `no`: still in 0.93-0;  X: disabled for 0.93-1


[^1]: Error in nls(formula = y ~ Expo(x, a, b), data = d.exp.Hlev, start = c(a = 1,  :
    step factor 0.000488281 reduced below 'minFactor' of 0.000976562

[^2]: Component "b1": Mean relative difference: 5.567653e-15
