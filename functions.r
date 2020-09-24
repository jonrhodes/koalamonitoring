smart.round <- function(x, digits = 0) {
# function to round while preserving the sum
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

projPops <- function(pops, nYears, decline, genYears, demogV, envtlV, envtlCorr, nSims) {
# function to build population simulations
# pops is a vector of starting population sizes for each population at t = 0
# nYears is how many years to run the simulation for
# decline in the simulated percent decline over 3 generations
# genYears is the number of years for three generations
# demogV is demographic stochasticity (variance in the growth rate among individuals)
# envtlV is the variance for environmental stochasticity
# envtlCorr is the correlation in environmental stochasticiity between populations
# nSims is how many independent simulations to run

  sim <- function(pops, nYears, r, demogV = 0, envtlVC = 0) {
  # function to run one iteration of the simulation
    Output <- round(matrix(pops))

    for (t in 1:nYears) {
      # draw random deviates and growth rate for this year
      # note here we adjust for the expectation of a log normal distribution: exp((sigma^2)/2)
      DemR <- rnorm(n = length(pops), mean = -(demogV / (Output[, ncol(Output)] ^ 2)) / 2, sd = sqrt(demogV) / pops)
      EnvR <- mvrnorm(n = 1, mu = -diag(envtlVC) / 2, Sigma = envtlVC)
      R <- rep(r, length(pops)) + DemR + EnvR

      Output <- cbind(Output, round(Output[, ncol(Output)] * exp(R)))
    }

    return(Output)
  }

  # construct the variance covariance matrix for environmental stochasiticity
  envtlVC <- diag(envtlV, nrow = length(pops), ncol = length(pops))
  envtlVC[lower.tri(envtlVC)] <- envtlV * envtlCorr
  envtlVC[upper.tri(envtlVC)] <- envtlV * envtlCorr

  # get r
  r <- ((1 - (decline / 100)) ^ (1 / genYears)) - 1

  # replicate simulations nSims times and return list
  return(replicate(n = nSims, expr = sim(pops = pops, nYears = nYears, r = r, demogV = demogV, envtlVC = envtlVC), simplify = FALSE))
}

getSurveys <- function(pops, survParams) {
# function to generate a survey sites and years
# pops is the population size data and the strata
# survParams is the survey parameters consisting of: f0, f0se, tranDayLine, tranDayArea, siteSize,
# budget, survIntens, method, monInter, strat, and monitRep

  # get area per site in ha
  areaSite <- as.numeric(survParams["survIntens"]) * as.numeric(survParams["siteSize"]) / 100
  # get line length per site in metres
  lengthSite <- (as.numeric(survParams["survIntens"]) * as.numeric(survParams["siteSize"]) * 100) / (2 / as.numeric(survParams["f0"]))
  # get days needed for each site
  daysPerSiteArea <- areaSite / as.numeric(survParams["tranDayArea"])
  daysPerSiteLine <- lengthSite / as.numeric(survParams["tranDayLine"])
  # get maximum number of sites per stratum
  pops_unn <- unnest(pops, cols = c(data)) %>% ungroup() %>% mutate(MAXSITES = AREA / as.numeric(survParams["siteSize"]))
  #get sums of weights for the different stratification strategies
  numStrat <- nrow(pops_unn)
  numStratUrb <- nrow(pops_unn[which(pops_unn$LUCategory == "Urban Footprint"),])
  numStratNonUrb <- numStrat - numStratUrb
  areaAll <- sum(pops_unn[, "AREA"])
  areaUrb <- sum(pops_unn[which(pops_unn$LUCategory == "Urban Footprint"), "AREA"])
  areaNonUrb <- areaAll - areaUrb
  DensAll <- sum(pops_unn[, "KD_MEAN"])
  DensUrb <- sum(pops_unn[which(pops_unn$LUCategory == "Urban Footprint"), "KD_MEAN"])
  DensNonUrb <- DensAll - DensUrb
  InvDensAll <- sum(1 / pops_unn[, "KD_MEAN"])
  InvDensUrb <- sum(1 / pops_unn[which(pops_unn$LUCategory == "Urban Footprint"), "KD_MEAN"])
  InvDensNonUrb <- InvDensAll - InvDensUrb

  # allocate sampling effort per stratum based on stratification procedure
  if (survParams["strat"] == "equal") {
    if (survParams["method"] == "line") {
      pops_unn <- pops_unn %>% mutate(SITES = as.numeric(survParams["budget"]) / (numStrat * daysPerSiteLine)) %>% mutate(SIZE = lengthSite, TYPE = "line")
      # re-allocate sites for strata with insufficient maximum number of sites
      if (length(which(pops_unn$SITES > pops_unn$MAXSITES)) > 0) {
        toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "MAXSITES"])
      } else {
        toAllocate <- 0
      }
      i <- 0 # note bails out if tries to allocate sites using more than 100 iterations
      while (toAllocate > 0 & i < 100) {
        allocSite <- toAllocate / length(which(pops_unn$SITES < pops_unn$MAXSITES))
        pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES), "SITES"] <- pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES), "SITES"] + allocSite
        pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "SITES"] <- pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "MAXSITES"]
        if (length(which(pops_unn$SITES > pops_unn$MAXSITES)) > 0) {
          toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "MAXSITES"])
        } else {
          toAllocate <- 0
        }
        i <- i + 1
        if (i == 100) {warning("more than 100 iterations to allocate survey effort")}
      }
      # round number of sites
      pops_unn$SITES <- smart.round(pops_unn$SITES)
    } else if ((survParams["method"] == "area")) {
      pops_unn <- pops_unn %>% mutate(SITES = as.numeric(survParams["budget"]) / (numStratNonUrb * daysPerSiteLine + numStratUrb * daysPerSiteArea)) %>% mutate(SIZE = lengthSite, TYPE = "line")
      pops_unn[which(pops_unn$LUCategory == "Urban Footprint"),"SIZE"] <- areaSite
      pops_unn[which(pops_unn$LUCategory == "Urban Footprint"),"TYPE"] <- "area"
      # re-allocate sites for strata with insufficient maximum number of sites (line transects)
      if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line")) > 0) {
        toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "MAXSITES"])
      } else {
        toAllocate <- 0
      }
      i <- 0 # note bails out if tries to allocate sites using more than 100 iterations
      while (toAllocate > 0 & i < 100) {
        allocSite <- toAllocate / length(which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "line"))
        pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] <- pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] + allocSite
        pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] <- pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "MAXSITES"]
        if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line")) > 0) {
          toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "MAXSITES"])
        } else {
          toAllocate <- 0
        }
        i <- i + 1
        if (i == 100) {warning("more than 100 iterations to allocate survey effort")}
      }
      # re-allocate sites for strata with insufficient maximum number of sites (area searches)
      if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area")) > 0) {
        toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "MAXSITES"])
      } else {
        toAllocate <- 0
      }
      i <- 0 # note bails out if tries to allocate sites using more than 100 iterations
      while (toAllocate > 0 & i < 100) {
        allocSite <- toAllocate / length(which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "area"))
        pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] <- pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] + allocSite
        pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] <- pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "MAXSITES"]
        if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area")) > 0) {
          toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "MAXSITES"])
        } else {
          toAllocate <- 0
        }
        i <- i + 1
        if (i == 100) {warning("more than 100 iterations to allocate survey effort")}
      }
      # round number of sites
      pops_unn[which(pops_unn$TYPE == "line"), "SITES"] <- smart.round(as.matrix(pops_unn[which(pops_unn$TYPE == "line"), "SITES"]))
      pops_unn[which(pops_unn$TYPE == "area"), "SITES"] <- smart.round(as.matrix(pops_unn[which(pops_unn$TYPE == "area"), "SITES"]))
    } else {
      stop("incorrect survey type")
    }
  } else if (survParams["strat"] == "area") {
    if (survParams["method"] == "line") {
      pops_unn <- pops_unn %>% mutate(SITES = AREA * as.numeric(survParams["budget"]) / (areaAll * daysPerSiteLine)) %>% mutate(SIZE = lengthSite, TYPE = "line")
      # re-allocate sites for strata with insufficient maximum number of sites
      if (length(which(pops_unn$SITES > pops_unn$MAXSITES)) > 0) {
        toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "MAXSITES"])
      } else {
        toAllocate <- 0
      }
      i <- 0 # note bails out if tries to allocate sites using more than 100 iterations
      while (toAllocate > 0 & i < 100) {
        allocSite <- toAllocate / length(which(pops_unn$SITES < pops_unn$MAXSITES))
        pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES), "SITES"] <- pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES), "SITES"] + allocSite
        pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "SITES"] <- pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "MAXSITES"]
        if (length(which(pops_unn$SITES > pops_unn$MAXSITES)) > 0) {
          toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "MAXSITES"])
        } else {
          toAllocate <- 0
        }
        i <- i + 1
        if (i == 100) {warning("more than 100 iterations to allocate survey effort")}
      }
      # round number of sites
      pops_unn$SITES <- smart.round(pops_unn$SITES)
    } else if ((survParams["method"] == "area")) {
      pops_unn <- pops_unn %>% mutate(SITES = AREA * as.numeric(survParams["budget"]) / (areaNonUrb * daysPerSiteLine + areaUrb * daysPerSiteArea)) %>% mutate(SIZE = lengthSite, TYPE = "line")
      pops_unn[which(pops_unn$LUCategory == "Urban Footprint"),"SIZE"] <- areaSite
      pops_unn[which(pops_unn$LUCategory == "Urban Footprint"),"TYPE"] <- "area"
      # re-allocate sites for strata with insufficient maximum number of sites (line transects)
      if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line")) > 0) {
        toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "MAXSITES"])
      } else {
        toAllocate <- 0
      }
      i <- 0 # note bails out if tries to allocate sites using more than 100 iterations
      while (toAllocate > 0 & i < 100) {
        allocSite <- toAllocate / length(which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "line"))
        pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] <- pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] + allocSite
        pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] <- pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "MAXSITES"]
        if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line")) > 0) {
          toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "MAXSITES"])
        } else {
          toAllocate <- 0
        }
        i <- i + 1
        if (i == 100) {warning("more than 100 iterations to allocate survey effort")}
      }
      # re-allocate sites for strata with insufficient maximum number of sites (area searches)
      if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area")) > 0) {
        toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "MAXSITES"])
      } else {
        toAllocate <- 0
      }
      i <- 0 # note bails out if tries to allocate sites using more than 100 iterations
      while (toAllocate > 0 & i < 100) {
        allocSite <- toAllocate / length(which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "area"))
        pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] <- pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] + allocSite
        pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] <- pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "MAXSITES"]
        if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area")) > 0) {
          toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "MAXSITES"])
        } else {
          toAllocate <- 0
        }
        i <- i + 1
        if (i == 100) {warning("more than 100 iterations to allocate survey effort")}
      }
      # round number of sites
      pops_unn[which(pops_unn$TYPE == "line"), "SITES"] <- smart.round(as.matrix(pops_unn[which(pops_unn$TYPE == "line"), "SITES"]))
      pops_unn[which(pops_unn$TYPE == "area"), "SITES"] <- smart.round(as.matrix(pops_unn[which(pops_unn$TYPE == "area"), "SITES"]))
    } else {
      stop("incorrect survey type")
    }
  } else if (survParams["strat"] == "density") {
    if (survParams["method"] == "line") {
      pops_unn <- pops_unn %>% mutate(SITES = KD_MEAN * as.numeric(survParams["budget"]) / (DensAll * daysPerSiteLine)) %>% mutate(SIZE = lengthSite, TYPE = "line")
      # re-allocate sites for strata with insufficient maximum number of sites
      if (length(which(pops_unn$SITES > pops_unn$MAXSITES)) > 0) {
        toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "MAXSITES"])
      } else {
        toAllocate <- 0
      }
      i <- 0 # note bails out if tries to allocate sites using more than 100 iterations
      while (toAllocate > 0 & i < 100) {
        allocSite <- toAllocate / length(which(pops_unn$SITES < pops_unn$MAXSITES))
        pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES), "SITES"] <- pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES), "SITES"] + allocSite
        pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "SITES"] <- pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "MAXSITES"]
        if (length(which(pops_unn$SITES > pops_unn$MAXSITES)) > 0) {
          toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "MAXSITES"])
        } else {
          toAllocate <- 0
        }
        i <- i + 1
        if (i == 100) {warning("more than 100 iterations to allocate survey effort")}
      }
      # round number of sites
      pops_unn$SITES <- smart.round(pops_unn$SITES)
    } else if ((survParams["method"] == "area")) {
      pops_unn <- pops_unn %>% mutate(SITES = KD_MEAN * as.numeric(survParams["budget"]) / (DensNonUrb * daysPerSiteLine + DensUrb * daysPerSiteArea)) %>% mutate(SIZE = lengthSite, TYPE = "line")
      pops_unn[which(pops_unn$LUCategory == "Urban Footprint"),"SIZE"] <- areaSite
      pops_unn[which(pops_unn$LUCategory == "Urban Footprint"),"TYPE"] <- "area"
      # re-allocate sites for strata with insufficient maximum number of sites (line transects)
      if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line")) > 0) {
        toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "MAXSITES"])
      } else {
        toAllocate <- 0
      }
      i <- 0 # note bails out if tries to allocate sites using more than 100 iterations
      while (toAllocate > 0 & i < 100) {
        allocSite <- toAllocate / length(which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "line"))
        pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] <- pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] + allocSite
        pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] <- pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "MAXSITES"]
        if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line")) > 0) {
          toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "MAXSITES"])
        } else {
          toAllocate <- 0
        }
        i <- i + 1
        if (i == 100) {warning("more than 100 iterations to allocate survey effort")}
      }
      # re-allocate sites for strata with insufficient maximum number of sites (area searches)
      if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area")) > 0) {
        toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "MAXSITES"])
      } else {
        toAllocate <- 0
      }
      i <- 0 # note bails out if tries to allocate sites using more than 100 iterations
      while (toAllocate > 0 & i < 100) {
        allocSite <- toAllocate / length(which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "area"))
        pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] <- pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] + allocSite
        pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] <- pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "MAXSITES"]
        if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area")) > 0) {
          toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "MAXSITES"])
        } else {
          toAllocate <- 0
        }
        i <- i + 1
        if (i == 100) {warning("more than 100 iterations to allocate survey effort")}
      }
      # round number of sites
      pops_unn[which(pops_unn$TYPE == "line"), "SITES"] <- smart.round(as.matrix(pops_unn[which(pops_unn$TYPE == "line"), "SITES"]))
      pops_unn[which(pops_unn$TYPE == "area"), "SITES"] <- smart.round(as.matrix(pops_unn[which(pops_unn$TYPE == "area"), "SITES"]))
    } else {
      stop("incorrect survey type")
    }
  } else if (survParams["strat"] == "invdensity") {
    if (survParams["method"] == "line") {
      pops_unn <- pops_unn %>% mutate(SITES = (1 / KD_MEAN) * as.numeric(survParams["budget"]) / (InvDensAll * daysPerSiteLine)) %>% mutate(SIZE = lengthSite, TYPE = "line")
      # re-allocate sites for strata with insufficient maximum number of sites
      if (length(which(pops_unn$SITES > pops_unn$MAXSITES)) > 0) {
        toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "MAXSITES"])
      } else {
        toAllocate <- 0
      }
      i <- 0 # note bails out if tries to allocate sites using more than 100 iterations
      while (toAllocate > 0 & i < 100) {
        allocSite <- toAllocate / length(which(pops_unn$SITES < pops_unn$MAXSITES))
        pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES), "SITES"] <- pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES), "SITES"] + allocSite
        pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "SITES"] <- pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "MAXSITES"]
        if (length(which(pops_unn$SITES > pops_unn$MAXSITES)) > 0) {
          toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES), "MAXSITES"])
        } else {
          toAllocate <- 0
        }
        i <- i + 1
        if (i == 100) {warning("more than 100 iterations to allocate survey effort")}
      }
      # round number of sites
      pops_unn$SITES <- smart.round(pops_unn$SITES)
    } else if ((survParams["method"] == "area")) {
      pops_unn <- pops_unn %>% mutate(SITES = (1 / KD_MEAN) * as.numeric(survParams["budget"]) / (InvDensNonUrb * daysPerSiteLine + InvDensUrb * daysPerSiteArea)) %>% mutate(SIZE = lengthSite, TYPE = "line")
      pops_unn[which(pops_unn$LUCategory == "Urban Footprint"),"SIZE"] <- areaSite
      pops_unn[which(pops_unn$LUCategory == "Urban Footprint"),"TYPE"] <- "area"
      # re-allocate sites for strata with insufficient maximum number of sites (line transects)
      if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line")) > 0) {
        toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "MAXSITES"])
      } else {
        toAllocate <- 0
      }
      i <- 0 # note bails out if tries to allocate sites using more than 100 iterations
      while (toAllocate > 0 & i < 100) {
        allocSite <- toAllocate / length(which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "line"))
        pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] <- pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] + allocSite
        pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] <- pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "MAXSITES"]
        if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line")) > 0) {
          toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "line"), "MAXSITES"])
        } else {
          toAllocate <- 0
        }
        i <- i + 1
        if (i == 100) {warning("more than 100 iterations to allocate survey effort")}
      }
      # re-allocate sites for strata with insufficient maximum number of sites (area searches)
      if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area")) > 0) {
        toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "MAXSITES"])
      } else {
        toAllocate <- 0
      }
      i <- 0 # note bails out if tries to allocate sites using more than 100 iterations
      while (toAllocate > 0 & i < 100) {
        allocSite <- toAllocate / length(which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "area"))
        pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] <- pops_unn[which(pops_unn$SITES < pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] + allocSite
        pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] <- pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "MAXSITES"]
        if (length(which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area")) > 0) {
          toAllocate <- sum(pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "SITES"] - pops_unn[which(pops_unn$SITES > pops_unn$MAXSITES & pops_unn$TYPE == "area"), "MAXSITES"])
        } else {
          toAllocate <- 0
        }
        i <- i + 1
        if (i == 100) {warning("more than 100 iterations to allocate survey effort")}
      }
      # round number of sites
      pops_unn[which(pops_unn$TYPE == "line"), "SITES"] <- smart.round(as.matrix(pops_unn[which(pops_unn$TYPE == "line"), "SITES"]))
      pops_unn[which(pops_unn$TYPE == "area"), "SITES"] <- smart.round(as.matrix(pops_unn[which(pops_unn$TYPE == "area"), "SITES"]))
    } else {
      stop("incorrect survey type")
    }
  } else {
    stop("incorrect stratification protocol")
  }

  # generate the survey sites and years in a #sites (M) * #time steps (T) matrix
  getSites <- function(x, monInter, monitRep) {
    maxInter <- ifelse(x["MAXSITES"] / x["SITES"] < 1, 1, floor(x["MAXSITES"] / x["SITES"]))[[1]]
    monInter <- min(monInter, maxInter)
    out <- matrix(NA, nrow = x["SITES"] * monInter, ncol = monitRep)
    for (i in 1:monitRep) {
      out[(((i %% monInter) * x["SITES"]) + 1):(((i %% monInter) * x["SITES"]) + x["SITES"]), i] <- i
    }
    return(out)
  }

  # generate survey
  pops_unn <- pops_unn %>% add_column(SURVEY = apply(dplyr::select(pops_unn, MAXSITES, SITES), MARGIN = 1, FUN = function (x, monInter, monitRep) {if(x["SITES"] > 0) {return(getSites(x, monInter, monitRep))} else {return(NULL)}},
                        monInter = as.numeric(survParams["monInter"]), monitRep = as.numeric(survParams["monitRep"])))

  return(pops_unn)
}
