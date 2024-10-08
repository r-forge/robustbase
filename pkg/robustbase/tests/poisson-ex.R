
library(robustbase)

## instead of relying on  system.file("test-tools-1.R", package="Matrix"):
doExtras <- robustbase:::doExtras()
SysI <- Sys.info()
## IGNORE_RDIFF_BEGIN
for(f in system.file("xtraR", c("test-tools.R",
                                "platform-sessionInfo.R"), # -> moreSessionInfo()
                     package = "robustbase", mustWork=TRUE)) {
    cat("source(",f,"):\n", sep="") ; source(f)
}
doExtras
moreSessionInfo(print. = "all")
(arch <- SysI[["machine"]]) # needed to distinguish platforms  further down
isMac <- SysI[["sysname"]] == "Darwin"
isSun <- SysI[["sysname"]] == "SunOS"
.M <- .Machine; str(.M[grep("^sizeof", names(.M))]) ## differentiate long-double..
noLD16 <- (.M$sizeof.longdouble != 16)
if(arch == "x86_64") {
    if(noLD16)
        arch <- paste0(arch, "--no-long-double")
    else if(osVersion == "Fedora 30 (Thirty)")
        arch <- paste0(arch, "_F30")
    # else keep  'arch'  unchanged
}
## IGNORE_RDIFF_END


#### Poisson examples from Eva Cantoni's paper

### Using Possum Data
### ================

data(possumDiv)

## Try to follow closely Cantoni & Ronchetti(2001), JASA
dim(X <- possum.mat[, -1]) # 151 13
str(y <- possum.mat[, "Diversity"])
##--- reduce the matrix from singularity ourselves:
X. <- possum.mat[, -c(1, match(c("E.nitens", "NW-NE"), colnames(possum.mat)))]
dim(X.)# 151 11

## "classical via robust: c = Inf :
Inf. <- 1e5 ## --- FIXME

## The following used to fail because glm.fit() returns NA coefficients
## now fine .. keep this as test!
glm.cr <- glmrob(y ~ X, family = "poisson", tcc = Inf.)
(scr <- summary(glm.cr))

scl <- summary(glm.cl <- glm   (Diversity ~ . , data=possumDiv, family=poisson))
sc2 <- summary(glm.c2 <- glmrob(Diversity ~ . , data=possumDiv, family=poisson, tcc = Inf.))
MMg <- model.matrix(glm.cl)

assert.EQ(coef(scl), coef(sc2), tol = 6e-6, giveRE=TRUE) # 1.37e-6
dnms <- list(colnames(MMg), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
cf.sc <- array(c(-0.9469439, 0.01192096, -0.2724059, 0.04022862, 0.03988606, 0.07173483,
                 0.01763833, -0.01534376, 0.1149216, 0.06675529, 0.1169463, -0.4889071,
                 ## SE
                 0.2655031, 0.02194661, 0.2859216, 0.01120463, 0.01438884, 0.03814053,
                 0.01059779, 0.1916126, 0.2724202, 0.1901612, 0.1902903, 0.2474653,
                 ## z val
                 -3.566603, 0.5431798, -0.9527294, 3.590356, 2.772014, 1.880803,
                 1.664341, -0.08007701, 0.421854, 0.3510457, 0.6145675, -1.975659,
                 ## P val
                 0.0003616393, 0.587006, 0.3407272, 0.0003302263, 0.00557107, 0.05999869,
                 0.09604432, 0.936176, 0.6731316, 0.7255541, 0.5388404, 0.04819339),
               dim = c(12L, 4L), dimnames = dnms)
assert.EQ(cf.sc, coef(sc2), tol = 4e-7, giveRE=TRUE) # 8.48e-8


## c = 2.0
summary(g2 <- glmrob(Diversity ~ . , data=possumDiv, family=poisson, tcc = 2.0, trace=TRUE))

## c = 1.6
glm.r <- glmrob(Diversity ~ . , data=possumDiv, family=poisson, tcc = 1.6, trace=TRUE)
(s.16 <- summary(glm.r))
str(glm.r)

## Now with *smaller* X (two variables less):
glm.c2 <- glmrob(y ~ X., family = "poisson", tcc = Inf.)
summary(glm.c2)

## c = 1.6,  x-weights, as in Cantoni-Ronchetti
glm.r2 <- glmrob(y ~ X., family = "poisson",
                 tcc = 1.6, weights.on.x = "hat")

## Now the same, for the direct possum data (no matrix),
## This indeed gives the same coefficients as in
## Table 3 of Cantoni+Ronchetti(2001): .. (tech.rep.):
glm.r2. <- glmrob(Diversity ~ ., family = "poisson", data=possumDiv,
                  tcc = 1.6, weights.on.x = "hat", acc = 1e-15)
## here iterate till convergence (acc = 10^(-15))

(sglm.r2 <- summary(glm.r2.))
## This is too accurate for S.E. (but we have converged to end)
cf2 <- matrix(c(-0.898213938628341, 0.269306882951903,
                0.00717220104127189, 0.0224349606070713,
                -0.25335520175528,  0.288588183720387,
                0.0403970350911325, 0.0113429514237665,
                0.0411096703375411, 0.0145996036305452,
                0.0730250489306713, 0.0386771060643486,
                0.0176994176433365, 0.0107414247342375,
                -0.0289935051669504,0.194215229266707,
                0.149521144883774,  0.271648514202971,
                0.0503262879663932, 0.191675979065398,
                0.0909870068741749, 0.192192515800464,
                -0.512247626309172, 0.250763990619973), 12,2, byrow=TRUE)
assert.EQ(cf2, unname(coef(sglm.r2)[, 1:2]), tol = 1e-9, giveRE=TRUE)#-> show : ~ 1.46e-11
stopifnot(abs(glm.r2.$iter - 18) <= 1) # 18 iterations on 32-bit (2008)

## MT estimator -- "quick" examples

if(!robustbase:::doExtras()) {
    cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
    quit()
}
## if ( doExtras ) -----------------------------------------------------

X1 <- cbind(1, X.)

if(FALSE) ## for debugging ...
    options(warn = 1, error=recover)
options(nwarnings = 1000) # def. 50

RNGversion("3.5.0") ## [TODO: adapt to "current" RNG settings]
set.seed(57)
showSys.time(
    ## m1 <- glmrobMT(x=X1, y=y)
    m1 <- glmrob(Diversity ~ ., data=possumDiv, family=poisson, method="MT")
)
summary(warnings())

stopifnot(m1$converged)
assert.EQ(m1$initial,
c(-0.851594294907422, -0.0107066895370536, -0.226958540075445, 0.0355906625338308,
  0.048010654640958, 0.0847493155436896, 0.0133604488401352, -0.024115201062159,
  0.0270535337324518, 0.146022135657894, -0.00751380783260833, -0.417638086169033)
          , tol = 1e-13, check.attributes=FALSE, giveRE=TRUE)


dput(signif(unname(coef(m1)), 11)) ## -->
## Something strange going on: R CMD check is different from interactive R, here.
## ???? [I see that the byte compiler is not listed in sessionInfo]
## In any case, take the dput(.) output from the *.Rout[.fail] file
## 2015-07-21: on 32-bit, the results *change* when re-run ???
##------- different results on different platforms

beta1 <- list(i686 =
## old florence:
## c(-0.83715700394, 0.0085488694315, -0.16734609346, 0.040965601691,
##   0.042387113444, 0.063146240793, 0.018632137866, -0.0062886781262,
##   0.11466679192, 0.091457894347, -0.025009954018, -0.66867971209)
## for a "moment": f32sfs-2; 2015-07-20
## c(-0.83818366695, 0.0085885492587, -0.1680548609, 0.040969491636,
##   0.042401438906, 0.063170238296, 0.018647880253, -0.0058039548495,
##   0.11500468542, 0.091940159895, -0.024804291737, -0.66861710581)
## f32sfs-2; 2015-07-21; in "R CMD check"/BATCH, *not* interactive
c(-0.83701057367, 0.0085408263511, -0.16692955779, 0.040980220489,
  0.042389760873, 0.063145608346, 0.018632314682, -0.0062819674369,
  0.11513144785, 0.091268054568, -0.025531439815, -0.66981350787)
## f32sfs-2; 2015-07-21, in R-devel, several times in a row:
## c(-0.83734949811, 0.008554484224, -0.16727333284, 0.040980350692,
##   0.042391751765, 0.06315585848, 0.018633222478, -0.0062978140762,
##   0.11509071086, 0.091463771235, -0.025113314023, -0.66955433495)
, "x86_64" =
c(-0.83723213945, 0.0085385261915, -0.16697112315, 0.040985126003,
  0.042400738973, 0.063168847366, 0.01863253681, -0.0064477807228,
  0.11488937188, 0.091283185006, -0.025627390293, -0.66995658693)
, "x86_64--no-long-double" =
c(-0.83710423989, 0.0085428949874, -0.16713845989, 0.040973904414,
  0.042391910971, 0.063159426394, 0.018629240073, -0.006362108938,
  0.1145563969, 0.091490891317, -0.025378427464, -0.66943593439)
, "x86_64_F30" =
c(-0.83703991366, 0.008536691385, -0.16707196217, 0.040980171987,
  0.042388781206, 0.063132162167, 0.018634264818, -0.0064298708197,
  0.11486525895, 0.091433901799, -0.025384338265, -0.66920847831)
)
## just FYI: difference 32-bit vs 64-bit:
assert.EQ(beta1[[1]], beta1[[2]], tol = 0.004, check.attributes=FALSE, giveRE=TRUE)
## Mean relative difference: 0.00142 [~ 2013-12]; 0.00273 [f32sfs-2; 2015-08]; then (R-devel 2015-07-21): 0.000916
assert.EQ(beta1[[2]], beta1[[3]], tol = 0.002, check.attributes=FALSE, giveRE=TRUE)
## Mean relative difference: 0.00082849  [2014-11]

## when bypassing BLAS in matprod()      vvvvv seen 0.001385 [Lx 64b]:
assert.EQ(coef(m1), beta1[[arch]], tol = 0.002, # typically 1e-10 is ok !!
          check.attributes=FALSE, giveRE=TRUE)

## The same, with another seed:
set.seed(64)
showSys.time(
    ## m2 <- glmrobMT(x=X1, y=y)
    m2 <- glmrob(Diversity ~ ., data=possumDiv, family=poisson, method="MT")
)
summary(warnings())

stopifnot(m2$converged)
if(FALSE)
dput(signif(unname(m2$initial), 13)) ## -->
assert.EQ(m2$initial, ## so this is *not* platform (32bit/64bit) dependent:
c(-1.204304813829, 0.02776038445201, -0.3680174045842, 0.04325746912892,
  0.03895315289169, 0.04537145479989, 0.02847987541025, 0.07073207523212,
  0.355491639539, 0.1822955449528, 0.1323720331562, -0.3419939094877)
          , tol = 1e-12, check.attributes=FALSE, giveRE=TRUE)

dput(signif(unname(coef(m2)), 11)) ## -->
beta2 <- list(i686 =
## florence?, Nov. 2014 (or even Dec 2013)
## c(-0.83698669149, 0.0085587296184, -0.16778044558, 0.040960021262,
##   0.042402954975, 0.063188868629, 0.018630275088, -0.0061015509403,
##   0.11385896307, 0.090966386294, -0.02572887737, -0.66945784056)
## f32sfs-2, July 2015, "R CMD .." (non-interactive!):
c(-0.83644647378, 0.0085365454367, -0.16770422458, 0.040958113098,
  0.04238796628, 0.063174324485, 0.018618360015, -0.0062357940483,
  0.11380146782, 0.090988141307, -0.025500338638, -0.66949122367)
## f32sfs-2, July 2015, interactive
## c(-0.83675287265, 0.0085383816807, -0.16763418359, 0.040968861778,
##   0.042399340988, 0.063148815999, 0.018624181637, -0.0061320761338,
##   0.11423331389, 0.0912474233, -0.025508101291, -0.66971416165)
, "x86_64" =
c(-0.83687097624, 0.0085341676033, -0.1674299545, 0.040968820903,
  0.042397459287, 0.063159075944, 0.018625582804, -0.0063140636571,
  0.11426134017, 0.091317308575, -0.025373078819, -0.66957444238)
, "x86_64--no-long-double" = # (2024-09: updated from Lx 64b-noLD)
c(-0.83708963633, 0.0085786151508, -0.16814236036, 0.040958368946,
  0.042399661659, 0.063166875696, 0.018634273269, -0.0056619915417,
  0.1142614614, 0.091055956969, -0.025688991294, -0.66932060023)
, "x86_64_F30" = ## Fedora 30, R-devel (2019-06-13):
c(-0.83651130836, 0.0085272636623, -0.16777225909, 0.040958046751,
  0.042398611622, 0.063169934556, 0.018622060538, -0.0067041556052,
  0.11358762483, 0.090950270043, -0.025393966426, -0.66916946118)
)
## just FYI: difference 32-bit vs 64-bit:
assert.EQ(beta2[[1]], beta2[[2]], tol = 0.001, check.attributes=FALSE, giveRE=TRUE)
## Mean relative difference: 0.0009487 [~2013-12]; .0008741591 [2024-09]
assert.EQ(beta2[[2]], beta2[[3]], tol = 0.004, check.attributes=FALSE, giveRE=TRUE)
## Mean relative difference: 0.0005119 [2014-11]; 0.0011933 2024-09

## when bypassing BLAS in matprod()      vvvvv seen 0.0002766 [Lx 64b], 0.001328 [Lx 64b-noLD]
assert.EQ(coef(m2), beta2[[arch]], tol = 0.0008 * (if(noLD16) 10 else 1), # typically 1e-10 is ok !!
          check.attributes=FALSE, giveRE=TRUE)
## slight changes of algorithm often change the above by ~ 4e-4 !!!

summary(warnings())

###---- Model Selection -----

## (not yet)  [ MM had this in ../../robGLM1/tests/quasi-possum.R ]

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
