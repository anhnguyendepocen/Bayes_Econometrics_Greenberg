# Exercise 3.10
rm(list = ls())
upLimN <- 500
n1 <- 50
n2 <- 50
m2 <- 10    # recaptures; m2 <= min(n1, n2)
nc <- n1 + n2    # n1 + n2 (total captures)
n <- nc - m2  # unique individuals captured
N <- max(1,n):upLimN
a <- 5
b <- 10
logPN <- lfactorial(N-1) - log(N) - lfactorial(N-n) + lgamma(2 * N - nc + b) - lgamma(2 * N + a + b)
PN1 <- exp(logPN)
sumPN <- sum(PN1)
PN <- PN1/sumPN
plot(N, PN, type = "h", ylab = "p(N)", bty = "l")