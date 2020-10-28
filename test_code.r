Test<-readRDS("output/popsims/popsim_decline0_envtlV0.01_envtlCorr0.rds")

mean(unlist(lapply(Test, FUN = function(x) {mean(x$V21 / x$V1)})))*100









Test<-fitStrataModel(data[[1]][[1]], pops_unn)

Test <- fitModels(data[[1]][1:3], pops)

Test3<-aggFits(Test$data[[1]])








out4[[i]] <- distsampOpen(lambdaformula = ~1, gammaformula = ~1, omegaformula = ~1, pformula = ~1, data = umf, K = max(out3[[i]], na.rm = T) * 2 + 1, keyfun = "halfnorm", output = "density", unitsOut = "ha", dynamics = "trend", method = "Nelder-Mead")

out4[[i]] <- pcountOpen(lambdaformula = ~1, gammaformula = ~1, omegaformula = ~1, pformula = ~1, data = umf, mixture = "P", K = max(out3[[i]], na.rm = T) * 2 + 1, dynamics = "trend", method = "Nelder-Mead")


Test <- distsampOpen(lambdaformula = ~1, gammaformula = ~1, omegaformula = ~1, pformula = ~1, data = data[[1]][[1]]$MODDAT[[1]], K = max(getY(data[[1]][[1]]$MODDAT[[1]]), na.rm = T)  + 1, keyfun = "halfnorm", output = "density", unitsOut = "ha", dynamics = "trend", method = "Nelder-Mead")

Test2 <- pcountOpen(lambdaformula = ~1, gammaformula = ~1, omegaformula = ~1, pformula = ~1, data = data[[5]][[1]]$MODDAT[[2]], mixture = "P", K = max(slot(data[[5]][[1]]$MODDAT[[2]],"y"), na.rm = T) * 2 + 1, dynamics = "trend", method = "Nelder-Mead")


class(data[[5]][[1]]$MODDAT[[1]])[1]



#Generate some data
  set.seed(123)
  gamma=0.8
  M=100; T=10; J=4 # M = sites, T = time steps, J = distance intervals
  L <- 4600
  dens <- 0.5
  f0 <- 0.0285

  sigma <- sqrt(2) / (f0 * sqrt(pi))

  y1 <- array(NA, c(M, J, T))
  N <- matrix(NA, M, T)
  S <- G <- matrix(NA, M, T-1)
  db <- c(0, 25, 50, 75, 100)

  #Half-normal, line transect
  g <- function(x, sig) {(sqrt(2) / (sig * sqrt(pi))) * exp(-x^2/(2*sig^2))}

  # get probabilities
  cp <- numeric(J)
  cp[1] <- integrate(g, db[1], db[2], sig=sigma)$value
  for(j in 2:J) {
    cp[j] <- integrate(g, db[j], db[j+1], sig=sigma)$value
  }
  cp[j+1] <- 1-sum(cp)

  for(i in 1:M) {
    N[i,1] <- rpois(1, (2 * L * dens) / (10000 * f0))
    y1[i,1:J,1] <- rmultinom(1, N[i,1], cp)[1:J]

    for(t in 1:(T-1)) {
        N[i,t+1] <- rpois(1, N[i,t]*gamma)
        y1[i,1:J,t+1] <- rmultinom(1, N[i,t+1], cp)[1:J]
        }
  }
  y <- matrix(y1, M)

  # Primary periods of observations
  # In this case there are no gaps
  primPer <- matrix(as.integer(rep(1:T, M)), nrow=M, ncol=T, byrow=TRUE)

  #Make a covariate
  sc <- data.frame(x1 = rnorm(M))

  umf <- unmarkedFrameDSO(y = y, siteCovs=sc, numPrimary=T, primaryPeriod=primPer, dist.breaks=db,
                          survey="line", unitsIn="m", tlength=rep(1, M))

  fit <- distsampOpen(lambdaformula = ~x1, gammaformula = ~1, omegaformula = ~1, pformula = ~1, data = umf, K = 50, keyfun="halfnorm", dynamics = "trend")

  # Empirical Bayes estimates of abundance for each site / year
  re <- ranef(fit)
  plot(re, layout=c(10,5), xlim=c(-1, 10))

=============================================

# Fake data
M <- 4 # number of sites
J <- 3 # number of distance classes
T <- 3 # number of primary periods

db <- c(0, 10, 20, 30) # distance break points

y <- matrix(c(
   5,4,3, 6,2,1, NA,NA,NA, # In bin 1: 5 detections in primary period 1, 6 in period 2
   0,0,0, 0,1,0, 2,1,3,
   2,1,1, 0,0,0, 2,2,4,
   1,1,0, 1,1,1, 0,1,0), nrow=M, ncol=J*T, byrow=TRUE)
y

# Primary periods of observations
# In this case there are no gaps
primPer <- matrix(as.integer(c(
    1,2,3,
    1,2,3,
    1,2,3,
    1,2,3)), nrow=M, ncol=T, byrow=TRUE)

#Site covs: M rows and 1 column per covariate
site.covs <- data.frame(x1=1:4, x2=factor(c('A','B','A','B')))
site.covs

#Yearly site covs on gamma/omega
ysc <- list(
    x3 = matrix(c(
        1,2,3,
        1,2,3,
        1,2,3,
        1,2,3), nrow=M, ncol=T, byrow=TRUE))

umf <- unmarkedFrameDSO(y=y, siteCovs=site.covs, yearlySiteCovs=ysc,
                        numPrimary=T, primaryPeriod=primPer,
                        dist.breaks=db, tlength = rep(400, M), survey="line", unitsIn="m")

umf                         # look at data
summary(umf)                # summarize

======================================================

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
