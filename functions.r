smart.round <- function(x, digits = 0) {
# function to round while preserving the sum
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

projPops <- function(pops, names, nYears, decline, genYears, demogV, envtlV, envtlCorr, nSims) {
# function to build population simulations
# pops is a vector of starting population sizes for each population at t = 0
# names is a vector of the names of the populations
# nYears is how many years to run the simulation for
# decline in the simulated percent decline over 3 generations
# genYears is the number of years for three generations
# demogV is demographic stochasticity (variance in the growth rate among individuals)
# envtlV is the variance for environmental stochasticity
# envtlCorr is the correlation in environmental stochasticiity between populations
# nSims is how many independent simulations to run

  sim <- function(pops, names, nYears, r, demogV, envtlVC) {
  # function to run one iteration of the simulation
    Output <- round(matrix(pops))

    for (t in 1:nYears) {
      # draw random deviates and growth rate for this year
      # note here we adjust for the expectation of a log normal distribution: exp((sigma^2)/2)
      DemR <- rnorm(n = length(pops), mean = -(demogV / (Output[, ncol(Output)] ^ 2)) / 2, sd = sqrt(demogV) / Output[, ncol(Output)])
      EnvR <- mvrnorm(n = 1, mu = -diag(envtlVC) / 2, Sigma = envtlVC)
      R <- rep(r, length(pops)) + DemR + EnvR

      Output <- cbind(Output, round(Output[, ncol(Output)] * exp(R)))
    }

    out <- as_tibble(cbind(as_tibble(names), as_tibble(Output))) %>% rename(POP = value)
    return(out)
  }

  # construct the variance covariance matrix for environmental stochasiticity
  envtlVC <- diag(envtlV, nrow = length(pops), ncol = length(pops))
  envtlVC[lower.tri(envtlVC)] <- envtlV * envtlCorr
  envtlVC[upper.tri(envtlVC)] <- envtlV * envtlCorr

  # get r
  r <- ((1 - (decline / 100)) ^ (1 / genYears)) - 1

  # replicate simulations nSims times and return list
  return(replicate(n = nSims, expr = sim(pops = pops, names = names, nYears = nYears, r = r, demogV = demogV, envtlVC = envtlVC), simplify = FALSE))
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

simData <- function(survey, popSim, f0, f0se, stripMissLow, stripMissMean, stripMissHigh) {
# function to generate survey data
# survey is the sites selection for a single monitoring parameter combination
# popSim is the replicates for a single population dynamics parameter combination
# f0 is the mean of the detection function
# f0se is the standard error of the detection function
# stripMissLow is the lower bound probability of missing a koala in a strip/area search
# stripMissMean is the mean probability of missing a koala in a strip/area search
# stripMissHigh is the upper bound probability of missing a koala in a strip/area search

  # calculate mean and se of f0 on the log scale (since we assume f0 is log-normally distributed)
  f0log <- log(f0) - (0.5 * (f0se^2))
  f0selog <- sqrt(log(exp(2 * log(f0) - log(f0se^2)) + 1))
  sim <- function(popSim, survey, f0, f0se, stripMissLow, stripMissMean, stripMissHigh) {
    # join simulations to survey for this iteration
    join <- left_join(x = survey, y = popSim, by = "POP")

    # output for densities
    out <- list()
    # output for counts
    out2 <- list()
    # output for perpendicular distances
    out3 <- list()
    # output for model data
    out4 <- list()

    # loop through populations
    for (i in 1:nrow(join)) {
      if (is.null(join$SURVEY[[i]])) {
        out[[i]] <- NA
        out2[[i]] <- NA
        out3[[i]] <- NA
        out4[[i]] <- NA
      } else {
        # get densities
        out[[i]] <- join$SURVEY[[i]]
        out[[i]][which(!is.na(out[[i]]))] <- as.numeric(join[i, paste("V", join$SURVEY[[i]][which(!is.na(out[[i]]))], sep="")]) / as.numeric(join[i,"V1"])
        out[[i]] <- out[[i]] * matrix(rep(as.numeric(join[i, "KD_MEAN"]), ncol(join$SURVEY[[i]])), nrow = nrow(join$SURVEY[[i]]), ncol = ncol(join$SURVEY[[i]])) # note here we assume densities in each site
                                                                                                                                                                 # within a stratum are the same - variation
                                                                                                                                                                 # among sites can be introduced via the standard
                                                                                                                                                                 # deviation (see below) but sometimes results in
                                                                                                                                                                 # unrealistically high densities
        # out[[i]] <- out[[i]] * matrix(rep(exp(rnorm(n = nrow(join$SURVEY[[i]]), mean = as.numeric(join[i, "KDLN_MEAN"]), sd = as.numeric(join[i, "KDLN_SD"]))), ncol(join$SURVEY[[i]])), nrow = nrow(join$SURVEY[[i]]), ncol = ncol(join$SURVEY[[i]]))
        # get counts and perpendicular distances
        if (join[i, "TYPE"] == "line") {
          f0Rand <- f0  # exp(rnorm(n = 1, mean = f0log, sd = f0selog)) note here we assume f0 is the same among strata and known with certainty, but uncertainty can be introduced via the standard error
          sigma <- sqrt(2) / (f0Rand * sqrt(pi))
          M <- nrow(out[[i]])
          T <- ncol(out[[i]])
          J <- 4
          db <- qhnorm(p = c(0, 0.25, 0.5, 0.75, 0.99), sigma = sigma)
          cp <- numeric(J)
          cp[1] <- phnorm(db[2], sigma = sigma)
          for(j in 2:J) {
            cp[j] <- phnorm(db[j + 1], sigma = sigma) - phnorm(db[j], sigma = sigma)
          }
          cp[j + 1] <- 1 - sum(cp)
          primPer <- matrix(as.integer(rep(1:T, M)), nrow = M, ncol = T, byrow = TRUE)
          # get counts
          out2[[i]] <- matrix(suppressWarnings(rpois(n = M * T, lambda = (2 * as.numeric(join[i, "SIZE"]) * out[[i]]) / (10000 * f0Rand))), nrow = M, ncol = T)
          # get perpendicular distances
          out3[[i]] <- array(NA, c(M, J, T))
          for (k in 1:M) {
            if (!is.na(out2[[i]][k, 1])) {
              out3[[i]][k, 1:J, 1] <- rmultinom(1, out2[[i]][k, 1], cp)[1:J]
            }
            for(t in 1:(T - 1)) {
                if (!is.na(out2[[i]][k, t + 1])) {
                  out3[[i]][k, 1:J, t + 1] <- rmultinom(1, out2[[i]][k, t + 1], cp)[1:J]
                }
            }
          }
          out3[[i]] <- matrix(out3[[i]], M)
          # get model data
          out4[[i]] <- unmarkedFrameDSO(y = out3[[i]], numPrimary = T, primaryPeriod = primPer, dist.breaks = db,
                                  survey = "line", unitsIn = "m", tlength = rep(as.numeric(join[i, "SIZE"]), M))
        } else if ((join[i, "TYPE"]) == "area") {
          missRand <- stripMissMean # runif(1, stripMissLow, stripMissHigh)  note here we assume detection error is is the same among strata and known with certainty,
                                      # but uncertainty can be introduced via the range of values
          M <- nrow(out[[i]])
          T <- ncol(out[[i]])
          J <- 1
          primPer <- matrix(as.integer(rep(1:T, M)), nrow = M, ncol = T, byrow = TRUE)
          # get counts
          out2[[i]] <- matrix(suppressWarnings(rbinom(n = M * T, size = round(as.numeric(join[i, "SIZE"]) * out[[i]]), prob = 1 - missRand)), nrow = M, ncol = T)
          # get perpendicular distances (not needed in this case)
          out3[[i]] <- NA
          # get model data
          out4[[i]] <- unmarkedFramePCO(y = out2[[i]], numPrimary = T, primaryPeriod = primPer)
        }
      }
    }

    out <- join %>% add_column(MODDAT = out4) %>% dplyr::select(ID, SIZE, MODDAT)

    return(out)
  }

  # get simulated data
  data <- lapply(popSim, FUN = sim, survey = survey, f0 = f0, f0se = f0se, stripMissLow = stripMissLow, stripMissMean = stripMissMean, stripMissHigh = stripMissHigh)

  return(data)
}

fitModels <- function(survData, pops) {
# survData is a list of survey data for all replicates for a given population scenario and given monitoring strategy

  pops_unn <- unnest(pops, cols = c(data)) %>% ungroup() %>% dplyr::select(POP, ID, AREA)

  #functions
  fitOneModel <- function (data) {
  # function to fit one model in one stratum

    if (!(class(data) == "unmarkedFrameDSO" & class(data) == "unmarkedFrameDSO")) {
      mod <- NA
      type <- NA
    } else {
      if (class(data) == "unmarkedFrameDSO") {
        mod <- distsampOpen(lambdaformula = ~1, gammaformula = ~1, omegaformula = ~1, pformula = ~1, data = data, K = max(getY(data), na.rm = T) * 2 + 5, keyfun = "halfnorm",
                output = "density", unitsOut = "ha", dynamics = "trend", method = "Nelder-Mead")
        type <- slot(mod, "fitType") #distsampOpen
      } else if (class(data) == "unmarkedFramePCO") {
        mod <- pcountOpen(lambdaformula = ~1, gammaformula = ~1, omegaformula = ~1, pformula = ~1, data = data, mixture = "P", K = max(getY(data), na.rm = T) * 2 + 5,
                dynamics = "trend", method = "Nelder-Mead")
        type <- slot(mod, "fitType") #"pcountOpen"
      } else {
        stop("incorrect survey type")
      }
    }

    # create output
    if (!(class(data) == "unmarkedFrameDSO" & class(data) == "unmarkedFrameDSO")) {
      out <- matrix(NA, nrow = 2, ncol = 3)
      dimnames(out)[[1]] <- c("mu", "sigma")
    } else {
      out <- rbind(coef(mod), SE(mod))
      dimnames(out)[[1]] <- c("mu", "sigma")
    }
    return(list(type, out))
  }

  aggFits <- function(fits) {
  # function to get weighted average trends estimates, standard errors and whether trend detected at p = 0.05 level
    extractTrend <- function(fit) {
      if (!is.na(fit[[1]])) {
        return(fit[[2]][1, 2])
      } else {
        return(NA)
      }
    }
    extractSE <- function(fit) {
      if (!is.na(fit[[1]])) {
        return(fit[[2]][2, 2])
      } else {
        return(NA)
      }
    }
    extractDens <- function(fit, size) {
      if (!is.na(fit[[1]])) {
        if (fit[[1]] == "distsampOpen") {
          return(exp(fit[[2]][1, 1] + (0.5 * fit[[2]][2, 1] ^ 2)))
        } else if (fit[[1]] == "pcountOpen") {
          return(exp(fit[[2]][1, 1] + (0.5 * fit[[2]][2, 1] ^ 2)) / size)
        } else {
          stop("wrong model type")
        }
      } else {
        return(NA)
      }
    }

    result <- fits %>% mutate(TREND = map_dbl(.x = FITS, .f = extractTrend), SE = map_dbl(.x = FITS, .f = extractSE), DENS = map2_dbl(.x = FITS, .y = SIZE, .f = extractDens)) %>%
                        mutate(DAREA = DENS * AREA)

    weights <- result$DAREA / sum(result$DAREA, na.rm = TRUE)

    out <- list()
    outtrend <-  sum(result$TREND * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
    outse <- sqrt(sum(result$SE ^ 2 * weights ^ 2, na.rm = TRUE) / (sum(weights, na.rm = TRUE) ^ 2))
    outp <- pnorm(outtrend, mean = 0, sd = outse)
    out$trend <- exp(outtrend + (0.5 * outse ^ 2))
    #out$se <- outse
    #out$p <- outp
    out$detect <- if (outp < 0.05) {1} else {0}

    return(out)
  }

  fitStrataModel <- function (stratData, pops_unn) {
  # function to fit a model to each stratum

    # fit models and return list
    fits <- lapply(stratData$MODDAT, FUN = fitOneModel)

    # add to data frame
    stratData <- stratData %>% dplyr::select(ID, SIZE) %>% add_column(FITS = fits)

    unnestData <- left_join(x = pops_unn, y = stratData, by = "ID")
    nestData <- unnestData %>% group_by(POP) %>% nest()

    # get whether a decline is detected for the region
    out1 <- lapply(list(unnestData), FUN = aggFits)
    out1 <- unlist(lapply(list(unnestData), FUN = aggFits))
    out1a <- out1[1] # trend
    out1b <- out1[2] # detect

    # get the number of declines detected for each population
    out2 <- lapply(nestData$data, FUN = aggFits)
    out2 <- unlist(out2)
    out2 <- out2[seq(2, length(out2), 2)] # select only the detections (not trends)
    out2 <- sum(out2)
    output <- c(out1a, out1b, out2)
    names(output) <- c("trend", "region", "allpops")

    return(output)
  }

  # fit models across all replicates
  listTest <- lapply(survData, FUN = fitStrataModel, pops_unn = pops_unn)


  # prepare output
  out <- matrix(unlist(listTest), nrow = 3, ncol = length(survData))
  dimnames(out) <- list(c("trend", "region", "pops"), NULL)
  trendSD <- sd(out["trend",])
  trendLow <- quantile(out["trend", ], 0.025)
  trendHigh <- quantile(out["trend", ], 0.975)
  out <- apply(out, MARGIN = 1, FUN = mean)
  out <- c(out, trendSD, trendLow, trendHigh)
  names(out) <- c("trend", "region", "pops", "trendsd", "trendlow", "trendhigh")

  return(out)
}
