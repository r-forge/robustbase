library(robustbase)

source(system.file("test-tools-1.R", package="Matrix", mustWork=TRUE))
## -> assertError(), showSys.time(), ...

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
stopifnot(all.equal(cfE, cf2, tol = 1e-9),
          abs(glm.r2.$iter - 18) <= 1) # 18 iterations on 32-bit (2008)

## MT estimator -- "quick" examples

if(!doExtras) {
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

stopifnot(m1$converged,
          all.equal(m1$initial,
 c(-0.82454076117497, -0.0107066895370536, -0.226958540075445, 0.0355906625338308,
   0.048010654640958, 0.0847493155436896, 0.0133604488401352,  -0.0270535337324515,
  -0.0511687347946107, 0.146022135657894, -0.00751380783260816, -0.417638086169032), tol = 1e-13)
          ,
          all.equal(m1$final,
 c(-0.722403147823884, 0.00852588186299624, -0.166744211963903, 0.0409725499046939,
   0.0423860322049144, 0.0631479521560575, 0.0186215385026053, -0.114400042988754,
   -0.120742808664024, 0.0914653880905431, -0.0253258042856295, -0.669220978140912), tol = 1e-13)
          ,
          TRUE)

## The same, with another seed:
set.seed(64)
showSys.time(
    m2 <- glmrobMT(x=X., y=y)
)

stopifnot(m2$converged,
          all.equal(m2$initial,
 c(-0.56423428799895, -0.00758749903943305, -0.0526308127315302, 0.0421793131527516,
   0.0332564905179959, 0.0665372292781289, 0.0136866325653848, -0.00801997988512953,
   -0.182495118468575, -0.00116493347439625, 0.0107950403793612, -0.730857299552309), tol = 1e-13)
,
          all.equal(m2$final,
 c(-0.722056443290325, 0.00854660241578943, -0.167247920449585, 0.0409774131757375,
   0.0423793067248048, 0.063126142186602,   0.0186346991795864, -0.114818413029325,
   -0.121245882858204, 0.091337074056828,  -0.0254800458522171, -0.669078421826111), tol = 1e-13)
,
          TRUE)


###---- Model Selection -----

## (not yet)  [ MM had this in ../../robGLM1/tests/quasi-possum.R ]

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
