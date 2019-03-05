set.seed(103)
G <- 500
draws <- 0
while (length(draws) < G + 1){
        U <- runif(3)
        if (U[2] <= exp(-log(U[1]) - ((-log(U[1]))^2)/2 - 1/2) & U[3] <= 1/2)
                 draws <- c(draws, -log(U[1]))
        if (U[2] <= exp(-log(U[1]) - ((-log(U[1]))^2)/2 - 1/2) & U[3] >  1/2)
                 draws <- c(draws,  log(U[1]))
 }
 target <- 2 + 2*draws[2:(G+1)]
 print(mean = ); mean(target))
 print(sd = ); sd(target))