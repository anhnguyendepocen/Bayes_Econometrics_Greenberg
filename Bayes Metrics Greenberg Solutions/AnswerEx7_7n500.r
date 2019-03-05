# source("AnswerEx7_7n500.r")

rm(list = ls())
graphics.off()
library(tmvtnorm)
library(coda)

# Create data.
set.seed(123)   # To get same values for output.

n <- 500

bet <- rgamma(1, 2, 1)
x <- rnorm(n, 0, 1)
u <- rnorm(n, 0, 1) 


y <- bet * x + u

sigxy <- sum(x * y)
sigxx <- sum(x * x)
bethat <- (sigxy - 1 + sqrt( (sigxy - 1)^2 + 4 * sigxx))/(2 * sigxx)

csf <- 1   # Tuning parameter for scale factor.
sf <- csf /(1/bethat^2 + sigxx)

# MH algorithm
G <- 5000
burnin <- 1000

betStore <- rep(0, burnin + G)
betStore[1] <- 1

betPost <- function(bet, y, x) {
	bet * exp(-bet) * exp(-(1/2) * sum((y - bet * x) * (y - bet * x)))
	}

	# Next line uses package tmvtnorm to generate burnin + G candidates.
bet.can <- 	rtmvt(burnin + G, mean = bethat, sigma = sf, df = 5, lower = 0, upper = Inf)
for (g in 2:(burnin + G)) {
	bet.cur <- betStore[g-1]
	postRatg <- betPost(bet.can[g], y, x)/betPost(bet.cur, y,  x)
	kernRatg <- dt(bet.cur, df = 5)/dt(bet.can[g], df = 5)
	alphaxyg <- postRatg * kernRatg 
	v <- runif(1)
	if (v < alphaxyg) betStore[g] <- bet.can[g] else betStore[g] <- betStore[g - 1]
	}
	
	# Compute acceptance rate.
acc <- 0	
betStore <- betStore[-(1:burnin)]	
betStore.mcmc <- as.mcmc(betStore)
for (g in 2:G) { if ( betStore[g] != betStore[g-1])  acc <- acc + 1 }
accRate <- format(acc/(G-1), digits = 3)

    # Outputs
cat(paste("MH acceptance rate = ", accRate, "\n"), sep = "")
cat(paste("True beta = ",format(bet, digits = 4), "\n"), sep = "")
cat(paste("Summary of MH run, n = ", n, "\n"))
print(summary(betStore.mcmc))

     # Graph

postscript(file = "C:/users/egreenberg/my documents/aamyfiles/egtex/egbookms2/AnswersToHomeworks/Ex7_7n500.eps", horizontal = F)	 
hist(betStore, 40, freq = F, main = paste("n = ", n), xlab = expression(beta), ylab = expression(paste(pi, "(", beta, "|y)")))
expression(paste(pi,"(",beta[3],"|",y,")"))
dev.off()