library(MCMCpack)
datain <- read.table("babiesI.data", header = T)

priorMeanBeta <- c(105, -10)  #  Assumes intercept is 105 ounces and 10 ounces less for smokers.
Var <- matrix(c(20, 0, 0, 10), nrow = 2, ncol = 2, byrow = T)
priorPrecBeta <- solve(Var)   # Prior precisions.
c0 <- 6   # Assumes mean of the error variance is 400 and variance of variance is 100.
d0 <- 100

# Note that MCMCpack automatically includes an intercept.
babies.out <- MCMCregress(bwt ~ smoke, b0 = priorMeanBeta, B0 = priorPrecBeta, c0 = c0, d0 = d0) 

summary(babies.out)
