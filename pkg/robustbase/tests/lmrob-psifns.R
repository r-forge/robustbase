### tests psi und tuning.psi, tuning.chi argument of lmrob.control

library(robustbase)

data(aircraft)

set.seed(1)
summary(mp0 <- lmrob(Y ~ ., data = aircraft, psi = 'bisquare', method = 'SMDM'))

set.seed(2)
summary(mp1 <- update(mp0, psi = 'optimal'))

set.seed(3)
summary(mp2 <- update(mp0, psi = 'ggw'))

set.seed(4)
summary(mp3 <- update(mp0, psi = 'welsh'))

set.seed(5)
summary(mp4 <- update(mp0, psi = 'ggw', tuning.psi = 5.85, tuning.chi = 5.5))

set.seed(6)
summary(mp5 <- update(mp0, psi = 'ggw', tuning.psi = 0.95, tuning.chi = 0.5))

