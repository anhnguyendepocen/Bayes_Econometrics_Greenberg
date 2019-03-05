# egRegression.r

   #  This is for regression with Gaussian errors.
egRegGauss <- function(y, X, b00, invB.0, sig2, a.0, d.0, burnIn, totIter)
{    #Begin function.
betamm <- matrix(0, nrow = totIter, ncol = dim(X)[2])
s2mm   <- rep(0, totIter)

# for g == 1
s2mm[1]    <- sig2
B1          <- solve(invB.0 + (1/s2mm[1])*t(X)%*%X)
meanBeta1   <- B1 %*% (b00 + (1/s2mm[1])*t(X)%*%y)
betamm[1,] <- rmvnorm(1, meanBeta1, B1)

for (g in 2:totIter)   {
	if ((g%%100)==0) cat("g = ", g, sep="", "\n")
	Bg <- solve(invB.0 + (1/s2mm[g-1])*t(X)%*%X)
	meanBetag <- Bg %*% (b00 + (1/s2mm[g-1])*t(X)%*%y)
	betamm[g,] <- rmvnorm(1, meanBetag, Bg)
	devg  <- y - X %*% betamm[g,]
	delg  <- (d.0 + t(devg) %*% devg)/2
	invs2 <- rgamma(1, (a.0 + n)/2, scale = 1/delg) 
	s2mm[g] <- 1/invs2           }

betam   <- betamm[(burnIn + 1):totIter, ]
s2m     <- s2mm[(burnIn + 1):totIter]   

return(list("betaStore" = betam, "s2Store" = s2m)) 
}

	# This is for regression with Student-t errors. 
	                     
egRegStudent <- function(y, X, b00, invB.0, alpha.1, d.0, nu1, burnIn, totIter)   {
betamm  <- matrix(0, nrow = totIter, ncol = dim(X)[2])
s2mm    <- rep(1, totIter)
lambdamm <- matrix(1, nrow = totIter, ncol = n)

# for g == 1
	#lambdamm[1, ] <- rep(1, n)
	Lambdamm1 <- diag(lambdamm[1, ])
	s2mm[1]    <- sig2
	B1       <- solve(invB.0 + (1/s2mm[1]) * t(X) %*% Lambdamm1 %*% X)
	meanBeta1  <- B1 %*% (b00 + (1/s2mm[1]) * t(X) %*% Lambdamm1 %*% y)
	betamm[1,] <- rmvnorm(1, meanBeta1, B1) 

for (g in 2:totIter) {
	if ((g %% 100)==0) cat("g = ", g, sep="", "\n")
# For beta.	   
	Lambdammg <- diag(lambdamm[g-1, ])
	Bg <- solve(invB.0 + (1/s2mm[g-1]) * t(X) %*% Lambdammg %*% X)
	meanBetag <- Bg %*% (b00 + (1/s2mm[g-1]) * t(X) %*% Lambdammg %*% y)
	betamm[g,] <- rmvnorm(1, meanBetag, Bg)
# For sigma^2	
	devg <- y - X %*% betamm[g, ]
	de1g <- (d.0 + t(devg) %*% Lambdammg %*% devg)/2
	invs2 <- rgamma(1, alpha.1/2, scale = 1/de1g) 
	s2mm[g] <- 1/invs2 
# For lambda	
	nu2g <- (nu0 + (1/s2mm[g]) * devg^2)/2
	lambdamm[g, ] <- rgamma(c(rep(1,n)), (nu1/2) * rep(1, n), scale = 1/nu2g)
}
betam   <- betamm[(burnIn + 1):totIter, ]
s2m     <- s2mm[(burnIn + 1):totIter]   
lambdam <- lambdamm[(burnIn + 1):totIter,]

return(list("betaStore" = betam, "s2Store" = s2m, "lambdaStore" = lambdam)) 
}

logMLt <- function(nu0, sig2Star, betaStar, y, X, b.0, B.0, alph.0, d.0, lambdaStore, sig2Store, totIter)
{
logLikelihood <- sum(log(gamma((nu0 +1 )/2))-log(sqrt(sig2Star)) - (1/2)*log(nu0 * pi) - log(gamma(nu0/2)) 
                       - ((nu0 + 1))/2 * log (1 + ((y-X %*% betaStar)^2)/(nu0 * sig2Star)))
                        
         # Priors		 		 
logPriorBeta <- dmvnorm(betaStar, mean = b.0, sigma = B.0, log = TRUE)
logPriorSig2 <- dgamma(1/sig2Star, alpha.0/2, scale = 2/d.0, log = TRUE)
invB.0 <- solve(B.0)
	# Denominator
		# Computation of pi(betaStar|y)
postBetaSum <- 0	
for (g in 1:G)  {
	if ((g %% 100)==0) cat("g = ", g, "  Calculating pi(betaStar|y)", sep="", "\n")
	Lambdag <- diag(lambdaStore[g, ])
	B1g <- solve((1/sig2Store[g]) * t(X) %*% Lambdag %*% X + invB.0)
	b1g <- B1g %*% ((1/sig2Store[g]) * t(X) %*% Lambdag %*% y + b.0)
	postBetag <- dmvnorm(betaStar, mean = b1g, sigma = B1g)
	postBetaSum <- postBetaSum + postBetag  }
logPostBeta <- log((1/G) * postBetaSum)

#		Reduced runs for postSig2Star.
dev <- y - X %*% betaStar	
dev2 <- dev^2
postSig2Sum <- 0 
# for g = 1
    s2Redmm  <- numeric(totIter)
	s2Redmm[1] <- sig2
	lambda1  <- rep(1, n)
	Lambda1  <- diag(lambda1)
	lambdamm <- matrix(rep(1,n), nrow = totIter, ncol = n)

for (g in 2:burnIn)  {
	if ((g %% 100)==0) cat("g = ", g, "  Reduced run, burnIn", sep="", "\n")
		# For lambda	
	nu2g <- (nu0 + (1/s2Redmm[g-1]) * dev2)/2
	lambdammg  <- rgamma(rep(1,n), (nu1/2) * rep(1,n), scale = 1/nu2g)
	lambdamm[g, ] <- lambdammg
	
		# For sigma^2
	Lambdammg <- diag(lambdamm[g,])
	delta1g	<- (d.0 + t(dev) %*% Lambdammg %*% dev)/2
	invs2 <- rgamma(1, alpha.1/2, scale = 1/delta1g) 
	s2Redmm[g] <- 1/invs2  }
	
for (g in (burnIn + 1):totIter) {
	if ((g %% 100)==0) cat("g = ", g, "  Reduced run", sep="", "\n")
	postSigSum <- 0 
		# For lambda	
	nu2g <- (nu0 + (1/s2Redmm[g-1]) * dev2)/2
	lambdammg  <- rgamma(rep(1,n), (nu1/2) * rep(1,n), scale = 1/nu2g)
	lambdamm[g, ] <- lambdammg
		# For sigma^2
	Lambdammg <- diag(lambdamm[g,])
	delta1g	<- (d.0 + t(dev) %*% Lambdammg %*% dev)/2
	invs2 <- rgamma(1, alpha.1/2, scale = 1/delta1g) 
	s2Redmm[g] <- 1/invs2  
	postSig2Sum  <- postSig2Sum + dgamma(1/sig2Star, alpha.1/2, scale = 1/delta1g) }

logPostSig2 <- log((1/G) * postSig2Sum)	
	
logMargLike <-  logLikelihood + logPriorBeta + logPriorSig2 - (logPostBeta + logPostSig2)  
return(logMargLike)
}