library(robust)

x <- rpois(100, lambda=1)

glmrob(x ~ 1, family=poisson())
