# This one does the t errors.

# source("AnswerEx8_3.r")
rm(list = ls())
graphics.off()
library(coda)
library(xtable)
library(mvtnorm)
source("egRegression.R")

dataInRaw <- read.table("babiesII.data", header = T)

# Remove incomplete observations.

noGest <- which(dataInRaw[,2] == 999)	
dataIn1 <- dataInRaw[-noGest, ]
noParity <- which(dataIn1[, 3] == 9)
# All observations have parity.
noAge <- which(dataIn1[, 4] == 99)
dataIn2 <- dataIn1[-noAge, ]
noHeight <- which(dataIn2[, 5] == 99)
dataIn3 <- dataIn2[-noHeight, ]
noWeight <- which(dataIn3[, 6] == 999)
dataIn4 <- dataIn3[-noWeight, ]
noSmoke <- which(dataIn4[, 7] == 9)
dataIn5 <- dataIn4[-noSmoke, ]

y <- dataIn5[, 1]
X <-as.matrix(dataIn5)
X1 <- cbind(X[, 2:7])
X2 <- cbind(X[, 2:6])

n <- length(y)
inter <- rep(1, n)
X1 <- cbind(inter, X1) 
X2 <- cbind(inter, X2)
k1 <- dim(X1)[2]
k2 <- dim(X2)[2]

  # Priors
b001 <- rep(0, k1)
b002 <- rep(0, k2)
B.01 <- diag(1, k1)
B.02 <- diag(1, k2)
invB.01 <- solve(B.01)
invB.02 <- solve(B.02)
alpha.0 <- 6
alpha.1 <- alpha.0 + n
d.0     <- 0.40   # Degrees of freedom for t distribution.
nu0    <- 5 
nu0vec <- nu0 * rep(1, n)
nu1    <- nu0 + 1
nu1vec <- nu1 * rep(1, n)
sig2 = 1   # Starting value for sigma2.

burnIn <- 1000
G <- 10000
totIter <- burnIn + G

   # Model 1t includes smoking varaiable.
model1t <- egRegStudent(y, X1, b001, invB.01, alpha.1, d.0, nu1, burnIn, totIter)  

betaStore1t <- model1t$betaStore
sig2Store1t <- model1t$s2Store    
lambdaStore1t <- model1t$lambdaStore  



# For marginal likelihood computations model 1
betaStar1t <- colMeans(betaStore1t)  
sig2Star1t <- mean(sig2Store1t)   

logMLmodel1t <- logMLt(nu0, sig2Star1t, betaStar1t, y, X1, b001, B.01, alpha.0, d.0, lambdaStore1t, sig2Store1t, totIter)


   # Model 2t excludes smoking variable.

model2t <- egRegStudent(y, X2, b002, invB.02, alpha.1, d.0, nu1, burnIn, totIter)

betaStore2t <- model2t$betaStore
sig2Store2t <- model2t$s2Store    
lambdaStore2t <- model2t$lambdaStore  

cat(paste("Summary, model with smoking, t errors"), "\n")
print(summary(cbind(betaStore1t, sig2Store1t)))

cat(paste("Summary, model without smoking, t errors"), "\n")
print(summary(cbind(betaStore2t, sig2Store2t)))

# For marginal likelihood computations model 2
betaStar2t <- colMeans(betaStore2t)  
sig2Star2t <- mean(sig2Store2t)  
lambdaStore2t <- model2t$lambdaStore  
logMLmodel2t <- logMLt(nu0, sig2Star2t, betaStar2t, y, X2, b002, B.02, alpha.0, d.0, lambdaStore2t, sig2Store2t, totIter)

cat(paste("logML, with smoking, t errors = ", logMLmodel1t, "\n"))
cat(paste("logML, without smoking, t errors", logMLmodel2t, "\n"))

dif.12t <- logMLmodel1t - logMLmodel2t
cat(paste("log10 BF12 = ", log10(exp(dif.12t))), "\n")

postscript(file = "C:/users/egreenberg/my documents/aamyfiles/egtex/egbookms2/AnswersToHomeworks/Ex8_3.eps", horizontal = F)
hist(betaStore1t[, 7], 40, freq = F, xlab = expression(beta[7]),  
	main =  expression(paste("Effect of smoking on birthweight, t errors: ", pi, "(", beta[7], "|y)")))
dev.off()
