rm(list = ls())
l <- 1:5
l <- l/10
s <- 1 - l
log10BF <- function(m, e, l, s) {
	n <- 2 * m 
	logN <- lgamma(e + 1) - (e + 1) * log(n  + 1)
	logD <- lgamma(e*l + 1) - (e*l + 1) * log(m + 1) + lgamma(e*s + 1) - (e*s + 1) * log(m + 1)
	logBF <- logN - logD
	log10BF <- log10(exp(logBF))
		}
BF100 <- log10BF(100, 10, l, s)
BF200 <- log10BF(200, 20, l, s)	