# source("AnswerEx7_8n500.r")

rm(list = ls())
graphics.off()
library(tmvtnorm)
library(coda)
library(pscl)  # This gives us the inverse gamma distritution.

# Create data.
set.seed(123)   # To get same values for output.

n <- 500

# The next two lines choose values of beta and sigma^2 from their prior distributions.
betTrue <- rgamma(1, 2, 1)
sig2True <- rigamma(1, 5/2, 3/2)  
x <- rnorm(n, 0, 1)
u <- rnorm(n, 0, 1) 


y <- betTrue * x + sqrt(sig2True) * u

sigxy <- sum(x * y)
sigxx <- sum(x * x)

csf <- 1   # Tuning parameter for scale factor.
sf <- function(csf, bet, sigma2)  csf /(1/bet^2 + sigxx/sigma2)

# MH algorithm
G <- 10000
burnin <- 1000

betStore <- rep(0, burnin + G)
sig2Store <- rep(0, burnin + G)
betStore[1] <- 1

betPost <- function(bet, sig2, y, x) {
	bet * exp(-bet) * exp(-(1/(2*sig2)) * sum((y - bet * x) * (y - bet * x)))
	}

a0 <- 5/2
b0 <- 3/2
a1 <- a0 + n/2
b1 <- function(bet, y, x) {3/2 + sum((y - bet * x)*(y - bet * x))/2  }

betaHat <- function(sig2) { num <- sigxy/sig2 -1 + sqrt((sigxy/sig2 - 1)^2 + 4 * sigxx/sig2)
                            den <- 2 * sigxx / sig2 
							return(num/den)   }

   # For g = 1, draw a sigma2.
b11 <- b1(betStore[1], y, x)
sig2Store[1] <- rigamma(1, a1, b11)

for (g in 2:(burnin + G)) {
      # First, beta given sigma^2
	bet.cur <- betStore[g-1]
	betMeang <- betaHat(sig2Store[g - 1])
	sfg <- sf(csf, betStore[g-1], sig2Store[g-1])
	betcang <- rtmvt(1, mean = betMeang, sigma = sfg, df = 5, lower = 0, upper = Inf)
	postRatg <- betPost(betcang, sig2Store[g-1], y, x)/betPost(bet.cur, sig2Store[g-1], y,  x)
	kernRatg <- dt(bet.cur, df = 5)/dt(betcang, df = 5)
	alphaxyg <- postRatg * kernRatg 
	v <- runif(1)
	if (v < alphaxyg) betStore[g] <- betcang else betStore[g] <- betStore[g - 1]
	  # Second, sigma^2 given beta
	b1g <- b1(betStore[g], y, x)
	sig2Store[g] <- rigamma(1, a1, b1g)
	}
sig2Store <- sig2Store[-(1:burnin)]
betStore  <- betStore[-(1:burnin)]	
	
	# Compute acceptance rate.
acc <- 0	

for (g in 2:G) { if ( betStore[g] != betStore[g-1])  acc <- acc + 1 }
accRate <- format(acc/(G-1), digits = 3)

    # Outputs
betStore.mcmc <- as.mcmc(betStore)
cat(paste("MH acceptance rate for beta = ", accRate, "\n"), sep = "")
cat(paste("True beta = ",format(betTrue, digits = 4), "\n"), sep = "")
cat(paste("Summary of MH run for beta, n = ", n, "\n"))
print(summary(betStore.mcmc))

sig2Store.mcmc <- as.mcmc(sig2Store)
cat(paste("True sigma^2 = ",format(sig2True, digits = 4), "\n"), sep = "")
cat(paste("Summary of run for sigma, n = ", n, "\n"))
print(summary(sig2Store.mcmc))

     # Graph

postscript(file = "C:/users/egreenberg/my documents/aamyfiles/egtex/egbookms2/AnswersToHomeworks/Ex7_8n500.eps", horizontal = F)	
op <- par()
op <-par(mfrow = c(2,2), pty = "s", cex.lab = 1.4)
	hist(betStore, 40, freq = F, main = paste("n = ", n), xlab = expression(beta), 
		ylab = expression(paste(pi, "(", beta, "|y)")))
	acf(betStore, ci = 0, main = expression(paste("AR for ", beta)))
	hist(sig2Store, 40, freq = F, main = paste("n = ", n), xlab = expression(sigma^2), 
		ylab = expression(paste(pi, "(", sigma^2, "|y)")))	
	acf(sig2Store, ci = 0, main = expression(paste("AR for ", sigma^2)))	
par(op)	
dev.off()