# G = the number of iterations.

Generator <- function(G=1000, burn.in=500){ 
    iter <- G + burn.in
    storage <- rep(NA, iter)
    # The initial state is set as State 1.
    x1 <- runif(1)
    if (x1  > 0.75) storage[1] <- 0 # 0 represents state 2
    if (x1 <= 0.75) storage[1] <- 1 # 1 represents state 1
    for (g in 2:iter){
    x <- runif(1)
    if (x >  0.125 & storage[(g-1)]==0) storage[g] <- 0
    if (x <= 0.125 & storage[(g-1)]==0) storage[g] <- 1
    if (x >  0.75  & storage[(g-1)]==1) storage[g] <- 0
    if (x <= 0.75  & storage[(g-1)]==1) storage[g] <- 1
    }
    return(storage[(burn.in+1):iter])
 }
 G  <- seq(from=2, to=5000, by=10)
 burn.in <- G # Here, burn.in is set to be equal to the iteration size.
 fraction <- rep(NA, length(G))
 for( i in 1:length(G)){
 fraction[i] <- sum(Generator(G=G[i], burn.in=burn.in[i]))/G[i]
 }