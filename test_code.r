#Generate some data
  set.seed(123)
  gamma=0.8
  M=100; T=10; J=4
  L <- 4600
  dens <- 0.05
  f0 <- 0.0285
  sigma <- sqrt(2) / (f0 * sqrt(pi))

  y <- array(NA, c(M, J, T))
  N <- matrix(NA, M, T)
  S <- G <- matrix(NA, M, T-1)
  db <- c(0, 25, 50, 75, 100)

  #Half-normal, line transect
  g <- function(x, sig) {(sqrt(2) / (sig * sqrt(pi))) * exp(-x^2/(2*sig^2))}

  # get probabilities
  cp <- numeric(J)
  cp[1] <- integrate(g, db[1], db[2], sig=sigma)$value
  for(j in 2:J) {
    a[j] <-  db[j+1]  - sum(a[1:j])
    cp[j] <- integrate(g, db[j], db[j+1], sig=sigma)$value
  }
  cp[j+1] <- 1-sum(cp)

  for(i in 1:M) {
    N[i,1] <- rpois(1, (2 * L * dens) / (10000 * f0))
    y[i,1:J,1] <- rmultinom(1, N[i,1], cp)[1:J]

    for(t in 1:(T-1)) {
        N[i,t+1] <- rpois(1, N[i,t]*gamma)
        y[i,1:J,t+1] <- rmultinom(1, N[i,t+1], cp)[1:J]
        }
  }
  y <- matrix(y, M)

  #Make a covariate
  sc <- data.frame(x1 = rnorm(M))

  umf <- unmarkedFrameDSO(y = y, siteCovs=sc, numPrimary=T, dist.breaks=db,
                          survey="line", unitsIn="m", tlength=rep(1, M))

  fit <- distsampOpen(lambdaformula = ~x1, gammaformula = ~1, omegaformula = ~1, pformula = ~1, data = umf, K = 50, keyfun="halfnorm", dymamics = "trend")

  # Empirical Bayes estimates of abundance for each site / year
  re <- ranef(fit)
  plot(re, layout=c(10,5), xlim=c(-1, 10))

  Test<-readRDS("popsim_decline30_envtlV0_envtlCorr0.5.rds")
  Data <- as.data.frame(t(Test[[1]]))
  names(Data) <- c("GoldCoast", "KoalaCoast", "MoretonBay", "Strad", "Noosa", "NWest", "SWest")

  ggplot() +
      geom_line(data = Data, aes(x = seq(0,50,1), y = GoldCoast), color = "red") +
      geom_line(data = Data, aes(x = seq(0,50,1), y = KoalaCoast), color = "blue") +
      geom_line(data = Data, aes(x = seq(0,50,1), y = MoretonBay), color = "green") +
      geom_line(data = Data, aes(x = seq(0,50,1), y = Strad), color = "purple") +
      geom_line(data = Data, aes(x = seq(0,50,1), y = Noosa), color = "black") +
      geom_line(data = Data, aes(x = seq(0,50,1), y = NWest), color = "yellow") +
      geom_line(data = Data, aes(x = seq(0,50,1), y = SWest), color = "orange") +
      xlab('Year') +
      ylab('Pop')

  Test<-readRDS("popsim_decline30_envtlV0.02_envtlCorr1.rds")
  Data <- as.data.frame(t(Test[[1]]))
  names(Data) <- c("GoldCoast", "KoalaCoast", "MoretonBay", "Strad", "Noosa", "NWest", "SWest")

      ggplot() +
          geom_line(data = Data, aes(x = seq(0,50,1), y = GoldCoast), color = "red") +
          geom_line(data = Data, aes(x = seq(0,50,1), y = KoalaCoast), color = "blue") +
          geom_line(data = Data, aes(x = seq(0,50,1), y = MoretonBay), color = "green") +
          geom_line(data = Data, aes(x = seq(0,50,1), y = Strad), color = "purple") +
          geom_line(data = Data, aes(x = seq(0,50,1), y = Noosa), color = "black") +
          geom_line(data = Data, aes(x = seq(0,50,1), y = NWest), color = "yellow") +
          geom_line(data = Data, aes(x = seq(0,50,1), y = SWest), color = "orange") +
          xlab('Year') +
          ylab('Pop')
