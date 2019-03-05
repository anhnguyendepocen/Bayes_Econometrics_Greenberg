set.seed(123)

rw <- function(p, q, T) { S <- rep(0, T)
	S[1] <- 0
	for (t in 2:T) {
		x <- runif(1)
		if (x <= p) S[t] <- S[t-1] + 1
		if (x > p & x <= (p + q)) S[t] <- S[t-1] - 1
		if (x > (p + q))  S[t] <- S[t-1] }
	return(S)	
	}
RW <- rw(0.3,0.4, 1000)

plot(RW, type = "l", xlab = "t", ylab = "State")