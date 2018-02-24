# Problems reported by CRAN (Brian Ripley) when checking robustbase with non-R versions of BLAS:

## Overview (manual, from the ~/R/Pkgs/robustbase/robustbase_<kind>.out  file):


| BLAS kind | file                  | error msg                  | error kind| comments |
|:--------- |:--------------------- |:----------                 |:----------| ---------|
| ATLAS     | tests/mc-strict.R     | a5$nonOut are not all TRUE | differing outlier set ||
| ATLAS     | tests/nlrob-tst.R     | assert.EQ(coef(nlR1),*, tol=1e-9) -- gave 1.073e-9 ||
| MKL       | tests/nlregrob-tst.R  | update(Cfit.40out, data = d.exp.Hlev) --> nls __error__ [^1] ||
| OpenBLAS  | tests/m-s-estimator.R | stopifnot(all.equal(mC[nm], mR[nm], ... tol=5e-15)) -> [^2] ||
| OpenBLAS  | tests/nlregrob-tst.R  | update(Cfit.40out, ==== same as MKL above .... ||


[^1]: Error in nls(formula = y ~ Expo(x, a, b), data = d.exp.Hlev, start = c(a = 1,  :
    step factor 0.000488281 reduced below 'minFactor' of 0.000976562

[^2]: Component "b1": Mean relative difference: 5.567653e-15
