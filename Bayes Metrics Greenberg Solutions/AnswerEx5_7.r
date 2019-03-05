set.seed(123)
G <- 1000
  M1  <- 1;  M2  <- -0.5
  S11 <- 2;  S12 <- 1
  S21 <- 1;  S22 <- 3
  A1  <-.2;  B1  <- 2
  A2  <- -1;  B2  <- 5
  a1    <- (A1 - M1)/sqrt(S11)
  b1    <- (B1 - M1)/sqrt(S11)

  # Sample Y1 using Inverse CDF method
  U     <-  runif(G)
  Y1    <-  qnorm(pnorm(a1) + U*(pnorm(b1) - pnorm(a1)))

  # Compute the Probability of interest using Importance Sampling and sampled Y1
  x1    <- Y1*sqrt(S11) + M1
  mu2   <- M2 + (S12/S11)*(x1 - M1)
  sigma22 <- S22 - S12^2/S11
  a2    <- (A2 - mu2)/sqrt(sigma22)
  b2    <- (B2 - mu2)/sqrt(sigma22)
  Pr    <- (1/G)*sum((pnorm(b2) - pnorm(a2))*(pnorm(b1) - pnorm(a1)))
  print(Pr)