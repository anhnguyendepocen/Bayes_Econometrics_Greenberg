# source("AnswerEx7_6.r")

rm(list = ls())
graphics.off()
set.seed(234)
alpha <- 5
beta  <- 5

G <- 1000
burnin <- 100
modeBeta <- (alpha - 1)/(alpha + beta -2)
c <- dbeta(modeBeta, alpha, beta) # 2.460938 (for alpha = beta = 5).
c <- c + .000001  # To be on the safe side.
oneOverc <- format(1/c, digits = 6)
storeAR <- 0.5
ctrAR <- 0
storeMH <- rep(0, G)
storeMH[1] <- 0.5

#  AR algorithm. 
while (length(storeAR) <= G) {
	ctrAR <- ctrAR + 1
	y <- runif(1)
	u <- runif(1)
	if (u <= dbeta(y, 5, 5)/c ) { storeAR <- c(storeAR, y)   }
	}

trueMean <- format(alpha/(alpha + beta), 3)
cat(paste("True Mean = alpha/(alpha + beta) = ", trueMean, "\n", "\n", sep = ""))	
		
meanAR <- format(mean(storeAR), 3)
cat(paste("Mean estimated by AR = ", meanAR, "\n", "\n", sep = ""))
accRateAR <- format(G/ctrAR, digits = 3)  # G divided by total number of trials needed to get G trials.	


#  MH algorithm.
for (g in 2:(G + burnin))
	{
	y <- runif(1)
	u <- runif(1)
	x <- storeMH[g-1]
	alphaxy <- dbeta(y, 5, 5)/dbeta(x, 5, 5)
	if (u <= alphaxy) {storeMH[g] <- y} else {storeMH[g] <- x}
	}
storeMH <- storeMH[-(1:burnin)]  # Eliminate burnin sample.
meanMH <- format(mean(storeMH), 3) 
cat(paste("Mean estimated by MH = ", meanMH, "\n", "\n", sep = ""))

ARHist <- hist(storeAR, 40)$intensities
MHHist <- hist(storeMH, 40)$intensities
maxY <- max(c(ARHist, MHHist))
yLimAR <- max(max(maxY), max(dbeta(x, 5, 5), from = 0, to = 1))
yLimMH <- max(max(maxY), max(dbeta(x, 5, 5), from = 0, to = 1))


	# To calculate acceptance rate for MH.
accMH <- 0
for (g in 2:G) { if (storeMH[g] != storeMH[g - 1]) accMH <- accMH + 1 } 
accRateMH <- format(accMH/G, 3)

cat(paste("AR acceptance rate = ", accRateAR, "  (1/c = ", oneOverc, ")",  "\n", "\n", sep = ""))
cat(paste("MH acceptance rate = ", accRateMH, "\n", "\n", sep = ""))

xx <- seq(from = 0, to = 1, length = 100)


postscript(file = "C:/users/egreenberg/my documents/aamyfiles/egtex/egbookms2/AnswersToHomeworks/Ex7_6.eps", horizontal = F)
op <- par()
op <- par(mfrow = c(2, 2), pty = "s", bty = "l", mar=c(2,4, 1 ,0.5), oma=c(1.5, 2, 1, 2))

# AR f(x) and cg(x)
curve(dbeta(x, 5, 5), main = "AR algorithm", ylim = c(0, yLimAR), ylab = "f(x), cg(x)")
lines(xx, rep(c, length(xx)), lty = 2)
legend(x = "topright", lty = c(1, 2), c("f(x)", "cg(x)"))

# AR hist of output and true.
hist(storeAR, 40, freq = FALSE, ylab = "f(x)", xlab = "x", main = "AR algorithm", ylim = c(0, yLimAR))
curve(dbeta(x, 5, 5), add = T)

curve(dbeta(x, 5, 5), main = "MH algorithm", ylim = c(0, yLimMH), ylab = "f(x), q(x)")
lines(xx, rep(1, length(xx)), lty = 2)
legend(x = "topright", lty = c(1, 2), c("f(x)", "q(x)"))

hist(storeMH, 40, freq = FALSE, ylab = "f(x)", xlab = "x", main = "MH algorithm", ylim = c(0, yLimMH))
curve(dbeta(x, 5, 5), add = T)		
par(op)	
dev.off()


