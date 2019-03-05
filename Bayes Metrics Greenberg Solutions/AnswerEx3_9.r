# Exercise 3.9
rm(list = ls())
upLimN <- 500
n <- 50
N <- (n+1):upLimN
a <- 5
b <- 10
logPN <- lfactorial(N-1) - log(N) - lfactorial(N-n) + lgamma(N - n + b) - lgamma(N + a + b)
PN1 <- exp(logPN)
sumPN <- sum(PN1)
PN <- PN1/sumPN
plot(N, PN, type = "h", ylab = "p(N)", bty = "l")