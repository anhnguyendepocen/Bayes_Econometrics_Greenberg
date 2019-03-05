rm(list = ls())
library(VGAM)	# For Laplace distribution.
library(coda)

set.seed(100)

Laplace_MH <- function(
 G,       # number of mcmc draws
 burnin,      # number of burnin draws
 x.start=0,       # initial value
 sigma)         # set the value of sigma
 {
 # initialize starting value and storage matrix
 x.holder <- matrix(NA, nrow = G + burnin, ncol = 1)
 x <- x.start
 # Metropolis-Hastings Loop
 for (g in 1:(G + burnin)) 
	{      # generate candidate
     x.can <- x + rnorm(1, mean=0, sd=sigma)
     # calculate acceptance probablity
     alpha <- min(exp(-abs(x.can) + abs(x)), 1)
     # accept or reject candidate based on alpha
     if (runif(1) < alpha) x <- x.can
     # store sample
     x.holder[g, 1] <- x
	}
 return(x.holder[-(1:burnin), 1])
 }
 
 # Print summaries.
 G <- 10000
 sigmav <- c(.05, 1, 2, 100)
 x.posterior <-matrix(0, nrow = G, ncol = length(sigmav))
 for (s in 1:4) { x.posterior[,s] <- Laplace_MH(G=10000, burnin=1000, sigma=sigmav[s]) 
                    print(paste("sigma = ", sigmav[s]));   print(summary(x.posterior[,s]));
   				}
lp25 <- format(qlaplace(.25), digits = 3); lp75 <- format(qlaplace(.75), digits = 3)				
print("True distribution")
paste("1st Qu. = ", lp25, "Median = Mean = 0 ", "3rd Qu. = ", lp75)

accRate <- numeric(4)
for (s in 1:length(sigmav)) {
	for(g in 2:G) {if (x.posterior[g, s] != x.posterior[g-1, s]) accRate[s] <- accRate[s] + 1 }
                    	 }
print("Acceptance Rate")
sigmavf <- format(sigmav, digits = 2)
print(sigmavf)
(accRate/G)

 # Make graphs.	
op <- par()
op <- par(mfrow = c(4, 2), pty = "s", bty = "l", mar=c(2,4, 1 ,0.5), oma=c(1.5, 2, 1, 2))
 for (s in 1:4){ zzz <- hist(x.posterior[,s], 20, plot = F)
                 zz <- max(zzz$count)   
                 hist(x.posterior[,s], 20, freq = F, xlim = c(-8, 8), xlab = "x", 
				 ylim = c(0, max(0.5, zz/G)), ylab = "f(x)", 
				 main = paste(expression(sigma), " = ", sigmav[s]) )
	             curve(dlaplace(x), -8, 8, add = T)			 
				 acf(x.posterior[,s])}
 par(op)	
 
 
