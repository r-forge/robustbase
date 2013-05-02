library(robustbase)

source(system.file("test-tools-1.R", package="Matrix", mustWork=TRUE))
## -> assertError(), showSys.time(), ...
## Newer version of the above test-tools-1.R contain this:
assert.EQ <- function(target, current, tol = if(show) 0 else 1e-15,
                      show = FALSE, ...) {
    ## Purpose: check equality *and* show non-equality
    ## ----------------------------------------------------------------------
    ## show: if TRUE, return (and hence typically print) all.equal(...)
    if(show) all.equal(target, current, tol = tol)
    else if(!isTRUE(r <- all.equal(target, current, tol = tol)))
	stop("all.equal() |-> ", paste(r, collapse=sprintf("%-19s","\n")))
}

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
## now fine
glm.cr <- glmrob(y ~ X, family = "poisson", tcc = Inf.)
(scr <- summary(glm.cr))

## c = 2.0
g2 <- glmrob(y ~ X, family = "poisson", tcc = 2.0, trace=TRUE)
summary(g2)

## c = 1.6
glm.r <- glmrob(y ~ X, family = "poisson", tcc = 1.6)
coef(summary(glm.r))

str(glm.r)

## Now with *smaller* X (two variable less):
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
cfE <- unname(coef(sglm.r2)[, 1:2])
all.equal(cfE, cf2, tol=0)#-> show : ~ 1.46e-11
assert.EQ(cfE, cf2, tol = 1e-9)
stopifnot(abs(glm.r2.$iter - 18) <= 1) # 18 iterations on 32-bit (2008)

## MT estimator -- "quick" examples

if(!robustbase:::doExtras()) {
    cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
    quit()
}
## if ( doExtras ) -----------------------------------------------------


if(FALSE) ## for debugging ...
options(warn = 1, error=recover)

set.seed(57)
showSys.time(
    m1 <- glmrobMT(x=X., y=y)
)

stopifnot(m1$converged)
assert.EQ(m1$initial,
 c(-0.82454076117497, -0.0107066895370536, -0.226958540075445, 0.0355906625338308,
   0.048010654640958, 0.0847493155436896, 0.0133604488401352,  -0.0270535337324515,
  -0.0511687347946107, 0.146022135657894, -0.00751380783260816, -0.417638086169032)
          , tol = 1e-13)

## MM: I'm shocked how much this changes after every tweak ...
beta1 <-
    c(-0.722624970605759, 0.00852749919992932, -0.166868774170829, 0.0409786388107873,
      0.0423939127370058, 0.0631485302047043, 0.0186282885623113, -0.11434877228984,
      -0.120681307761841, 0.091425088547055, -0.0253897328183308, -0.669546401877007)
beta1 <-
    c(-0.722420662687891, 0.0085476612701357, -0.166923595194219, 0.0409726077990677,
      0.0423823207249915, 0.0631502421922883,  0.018628347957457, -0.114433619812227,
      -0.120518637621447, 0.0912580486966675, -0.0256851944443132,-0.66933007259972)

all.equal(m1$final, beta1, tol = 0)#-> *see* the diff
assert.EQ(m1$final, beta1, tol = 1e-10)

## The same, with another seed:
set.seed(64)
showSys.time(
    m2 <- glmrobMT(x=X., y=y)
)

stopifnot(m2$converged)
assert.EQ(m2$initial,
c(-0.848813174290401, 0.0277603844520112, -0.368017404584204, 0.0432574691289175,
  0.0389531528916901, 0.0453714547998864, 0.0284798754102488, -0.355491639538963,
  -0.284759564306845, 0.182295544952766, 0.132372033156218, -0.34199390948774)
          , tol = 1e-13)

beta2 <-
 c(-0.722466858080995, 0.00853214741170862, -0.167080756780928, 0.0409708489859932,
   0.0423892852206197, 0.0631575676559838,  0.0186247550404213, -0.114368960856547,
   -0.12071086045638,  0.0913163131743179, -0.0254084789601603, -0.669355337605894)

beta2 <-
 c(-0.719639958181682, 0.00850921638739351, -0.166172980617748, 0.0410112563931713,
   0.0423638710960024, 0.0630362342605086, 0.0186359589765251, -0.116645159395719,
 -0.123130115061652, 0.0910865027225399, -0.0256449044169698, -0.67024227284216)

all.equal(m2$final, beta2, tol = 0)#-> *see* the diff; often in the order 4e-4
assert.EQ(m2$final, beta2, tol = 1e-10)

###---- Model Selection -----

## (not yet)  [ MM had this in ../../robGLM1/tests/quasi-possum.R ]

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
