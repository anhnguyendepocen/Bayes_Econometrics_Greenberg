# Francis Vella and Marno Verbeek, "Whose Wages Do Unions Raise? A
# Dynamic Model of Unionism and Wage Rate Determination for Young Men",
# Journal of Applied Econometrics, Vol. 13, No. 2, 1998, pp. 163-183.

# The data are taken from the National Longitudinal Survey (NLS Youth
# Sample) and contain observations on 545 males for the years 1980-1987.
# The variables are in alphabetical order, except for the first two which
# indicate the NLS individual identification number and the year of
# observation. The observations are ordered first by individual and
# second by period. The variables are saved in column space delimited
# format. Each line in the file contains one observation (NR, YEAR). 

# The data are in the DOS file vv-data.dat, which is zipped in
# vv-data.zip.

# See Table 1 of the paper for a description of the variables and summary
# statistics. 

# Variables:
# NR YEAR AG BLACK BUS CON ENT EXPER FIN HISP HLTH HOURS MAN MAR MIN NC NE 
# OCC1 OCC2 OCC3 OCC4 OCC5 OCC6 OCC7 OCC8 OCC9 PER PRO PUB RUR S SCHOOL TRA 
# TRAD UNION WAGE 

rm(list = ls())

library(MCMCpack)
library(xtable)

vvData   <- read.table("VV-DATA.DAT")

vvData87ind <- which(vvData[,2] == 1987)
vvData87    <- vvData[vvData87ind,]

y <- vvData87[, 36] # Defines dependent variable (log wage).
n <- length(y)

choosecov = c(3:11, 13:25, 27,28, 30:35)
X <- vvData87[,choosecov]
X <- as.matrix(X)
k <- dim(X)[2]

# Note that MCMCpack will supply an intercept.

# Note:  in output, variable 1 is intercept, variable 31 is union dummy, and variable 36 is sigma^2.
#        In summary, union is variable 35.

b.0      <- rep(0, k+1)    # Prior:  mean of beta in b0, variance in B0.
b.0P <- c(0, .10, .20)


B.0 <-  diag(nrow = k+1, ncol = k+1)
B.0P <- c(.0010, .0036, .0050)
# B.0[k+1, k+1] <- .0036
# B.0 <- solve(B.0)  # In MCMCpack, B.0 is precision matrix.

c.0 <- 6     # Parameters for sigma^2 in c0 and d0.
d.0 <- 0.40

outIJ <- matrix(0, nrow = length(b.0P), ncol = length(B.0P))
set.seed(100)
for (i in 1:3)   {  b.0[k+1] <-  b.0P[i]
	for (j in 1:3) { B.0[k+1, k+1] <- B.0P[j]; B.0I <- solve(B.0)
vv.reg.out <- MCMCregress(y ~ I(X), b0 = b.0, B0 = B.0I, c0 = c.0, d0 = d.0, 
    marginal.likelihood = "Chib95") 
	meanUnion <- format(mean(vv.reg.out[, 31]), digits = 3)
    cat(paste("Mean of beta7 = ", meanUnion, "Prior mean = ", b.0[k+1]), "Prior var =  ", B.0[k+1, k+1], "\n") 
			outIJ[i, j] <- meanUnion		}     }  	

  # To make a table of results.
  rownames(outIJ) <- c(b.0P)
  colnames(outIJ) <- c(B.0P)
  outIJ <- xtable(outIJ, caption = "Mean of $\\beta_U$", digits =3)
  print(outIJ)


