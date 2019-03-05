# source("AnswerEx8_1.r")

rm(list = ls())
graphics.off()
library(MCMCpack)
set.seed(123)

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

    # Model 1:  with smoking variable.

b01 <- c(105, 0, 0, 0, 0, 0, -10)
var1 <- diag(c(20, 10, 10, 10, 10, 10, 10))
B01 <- solve(var1)
c0 <- 6
d0 <- 100
model1 <- MCMCregress(y ~ X1, b0 = b01, B0 = B01, c0, d0, marginal.likelihood = "Chib95")
cat("Summary, model with smoking, Gaussian errors")
print(summary(model1))

   # Model 2:  without smoking variable.
   
b02 <- c(105, 0, 0, 0, 0, 0)
var2 <- diag(c(20, 10, 10, 10, 10, 10))
B02 <- solve(var2)
model2 <- MCMCregress(y ~ X2, b0 = b02, B0 = B02, c0, d0, marginal.likelihood = "Chib95")
cat("Summary, model without smoking, Gaussian errors")
print(summary(model2))

logmarglike1 <- attr(model1, "logmarglike")			  
logmarglike2 <- attr(model2, "logmarglike")			  
cat(paste("Normal errors, logMLModel 1 = ", logmarglike1, "\n"))
cat(paste("Normal errors, logMLModel 2 = ", logmarglike2, "\n"))

dif.12 <- logmarglike1 - logmarglike2
cat(paste("Normal errors, log10 BF12 = ", log10(exp(dif.12))), "\n")

# postscript(file = "C:/users/egreenberg/my documents/aamyfiles/egtex/egbookms2/AnswersToHomeworks/Ex8_1.eps", horizontal = F)
 hist(model1[, 7], 40, freq = F, xlab = expression(beta[7]),  
	main =  expression(paste("Effect of smoking on birthweight, Gaussian errors: ", pi, "(", beta[7], "|y)")))
#dev.off()