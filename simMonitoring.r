# load required libraries
library(tidyverse)
library(unmarked)
library(MASS)
library(parallel)
library(extraDistr)
library(binom)
library(abind)

# location of functions
source("functions.r")

# define input parameters
genYears <- 20 # number of years for a generation
demogV <- 0.4373 # demographic stochasticity variance - this is variance in individual annual growth rate among individuals - taken from Rhodes et al. (2011) based on the radio tacking data only
                # model. Here the standard deviation of the growth rate parameter estimate was 0.04458435 from 220 individuals, so we calculated demogV as (0.04458435^2) * 220
envtlV <- c(0, 0.01) # environmental stochasticity variance
envtlCorr <- c(0, 0.5, 1) # correlation in environmental stochasticity among populations
nSimsPop <- 1000 # number of replicates for the population simulations
nYears <- 50 # simulation time horizon
decline <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50) # percentage decline over three generations
f0 <- 0.0285 # distance detection parameter from Dique et al(2003), eff(W) = 1/f0
f0se <- 0.0021 # standard error of f0 from Dique et al(2003)
stripMissLow <- 0 # low end of range of proportion missed for strip transects from Dique et al. (2001)
stripMissMean <- 0.16 # mean proportion missed for strip transects from Dique et al. (2001)
stripMissHigh <- 0.33 # high end of range of proportion missed for strip transects from Dique et al. (2001)
tranDayLine <- 4600 # line transect length per day in m
tranDayArea <- 7.5 # area search area per day in ha
siteSize <- 100 # site size in ha
monitRep <- 21 # monitoring time frame
budget <- c(95, 190, 380) # field days per year
survIntens <- c(30, 60, 90) # proportion of the site surveyed
method <- c("line", "area") # "line" is all line and "area" is line transects outside the urban footprint and all area searches inside the urban footprint
monInter <- c(1, 2) # monitoring interval in years
strat <- c("equal", "area", "density", "invdensity")

# load strata and initial population densities - set -9999 to NA and remove strata with < 100 ha of habitat
# then group by genetic population and nest
# note here that we replace the density data for "core - urban footprint - north west" stratum with the density data for "core - urban footprint - north west" stratum
# since this density data is missing for "core - urban footprint - north west". we also use data from Biolink 2019 for North Stradbroke Island (NSI) and assume core and
# non-core habitat densities are the same for NSI. this is the file "start_dens_modified.csv"
pops <- read.csv("input/start_dens_modified.csv") %>%
  mutate(KD_MEAN = replace(KD_MEAN, KD_MEAN == -9999, NA), KD_SD = replace(KD_SD, KD_SD==-9999, NA),
    KDLN_MEAN = replace(KDLN_MEAN, KDLN_MEAN == -9999, NA), KDLN_SD=replace(KDLN_SD, KDLN_SD == -9999, NA)) %>%
  filter(AREA >= 100) %>% mutate(ABUND = round(KD_MEAN * AREA)) %>% as_tibble() %>% rowid_to_column("ID") %>% group_by(POP) %>% nest()
# calculate the population size for each genetic populations
pops <- pops %>% mutate(N = map_dbl(data, ~sum(.$ABUND)))

# simulate population dynamics - comment out if already run

procTime <- system.time({

popSimComb <- expand.grid(decline = decline, envtlV = envtlV, envtlCorr = envtlCorr)
# remove rows when environmental variation is 0 and correlation is 0 or 0.5
popSimComb <- popSimComb %>% filter(!(envtlV == 0 & ((envtlCorr == 0) | (envtlCorr == 0.5))))
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, list("pops", "projPops", "nSimsPop", "nYears", "genYears", "demogV" ))
clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(MASS))
parApply(cl = cl, X = popSimComb, MARGIN = 1, FUN = function(x) {saveRDS(object = projPops(pops = pops$N, names = pops$POP, nYears = nYears, decline = x["decline"], genYears = genYears,
                demogV = demogV, envtlV = x["envtlV"], envtlCorr = x["envtlCorr"], nSims = nSimsPop), file = paste("output/popsims/popsim", "_decline",
                x["decline"], "_envtlV", x["envtlV"], "_envtlCorr", x["envtlCorr"], ".rds", sep = ""))})
stopCluster(cl)

})

print(procTime)

# generate survey designs - comment out if already run

procTime <- system.time({

# get combinations of parameters
survSimComb <- expand.grid(budget = budget, survIntens = survIntens, method = method, monInter = monInter, strat = strat, monitRep = monitRep)
# collate fixed survey parameters together
fixedSurvParams <- c(f0, f0se, tranDayLine, tranDayArea, siteSize)
names(fixedSurvParams) <- c("f0", "f0se", "tranDayLine", "tranDayArea", "siteSize")
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, list("getSurveys", "smart.round", "pops", "fixedSurvParams" ))
clusterEvalQ(cl, library(tidyverse))
parApply(cl = cl, X = survSimComb, MARGIN = 1, FUN = function(x) {saveRDS(object = getSurveys(pops = pops, survParams = c(x, fixedSurvParams)), file = paste("output/surveys/survey", "_budget",
                x["budget"], "_survIntens", as.numeric(x["survIntens"]), "_method", x["method"], "_monInter", as.numeric(x["monInter"]), "_strat", x["strat"], "_monitRep", as.numeric(x["monitRep"]), ".rds", sep = ""))})
stopCluster(cl)

})

print(procTime)

# generate simulated survey data
# get survey information
Surveys <- apply(survSimComb, MARGIN = 1, FUN = function(x) {readRDS(file = paste("output/surveys/survey", "_budget", x["budget"], "_survIntens",
              as.numeric(x["survIntens"]), "_method", x["method"], "_monInter", as.numeric(x["monInter"]), "_strat", x["strat"], "_monitRep",
              as.numeric(x["monitRep"]), ".rds", sep = ""))})
# loop through population combinations
for (i in 1:nrow(popSimComb)) {

  procTime <- system.time({

  # get population simulations
  popSim <- readRDS(file = paste("output/popsims/popsim", "_decline", popSimComb[i, "decline"], "_envtlV",
        popSimComb[i, "envtlV"], "_envtlCorr", popSimComb[i, "envtlCorr"], ".rds", sep = ""))
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, list("simData", "popSim", "f0", "f0se", "stripMissLow", "stripMissMean", "stripMissHigh"))
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(extraDistr))
  clusterEvalQ(cl, library(unmarked))
  dataTemp <- parLapply(cl = cl, X = Surveys, fun = function(x) {return(simData(x, popSim, f0, f0se, stripMissLow, stripMissMean, stripMissHigh))})
  stopCluster(cl)
  saveRDS(object = dataTemp, file = paste("output/data/simdata", "_decline", popSimComb[i, "decline"], "_envtlV",
        popSimComb[i, "envtlV"], "_envtlCorr", popSimComb[i, "envtlCorr"], ".rds", sep = ""))
  rm(popSim)
  rm(dataTemp)
  gc()

  })

  print(procTime)

}

# fit models
# loop through population combinations
for (i in 1:nrow(popSimComb)) {

  procTime <- system.time({

  # get simulated data
  surveyData <- readRDS(file = paste("output/data/simdata", "_decline", popSimComb[i, "decline"], "_envtlV",
        popSimComb[i, "envtlV"], "_envtlCorr", popSimComb[i, "envtlCorr"], ".rds", sep = ""))
  # pick how many replicates - edit out if using all replicates
  surveyData <- lapply(surveyData, FUN = function(x) {x[1:100]})
  expR <-  log((1 - (popSimComb[i, "decline"] / 100)) ^ (1 / genYears))
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, list("fitModels", "pops", "expR"))
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(unmarked))
  clusterEvalQ(cl, library(binom))
  resultTemp <- parLapply(cl = cl, X = surveyData, fun = function(x) {return(fitModels(x, pops, expR))})
  stopCluster(cl)
  saveRDS(object = resultTemp, file = paste("output/results/result", "_decline", popSimComb[i, "decline"], "_envtlV",
                popSimComb[i, "envtlV"], "_envtlCorr", popSimComb[i, "envtlCorr"], ".rds", sep = ""))
  rm(surveyData)
  rm(resultTemp)
  gc()

  })

  print(procTime)

}

# compile results and output rankings

# add row ID to survey combinations
survSimComb <- survSimComb %>% rowid_to_column("ID")

# loop through population simulation combinations
for (i in 1:nrow(popSimComb)) {
  # get type 1 error
  type1 <- readRDS(file = paste("output/results/result", "_decline", 0, "_envtlV",
      popSimComb[i, "envtlV"], "_envtlCorr", popSimComb[i, "envtlCorr"], ".rds", sep = ""))
  # get results
  result <- readRDS(file = paste("output/results/result", "_decline", popSimComb[i, "decline"], "_envtlV",
      popSimComb[i, "envtlV"], "_envtlCorr", popSimComb[i, "envtlCorr"], ".rds", sep = ""))

  if (popSimComb[i, "decline"] > 0) {

    # compile type 1 errors
    type1Mat <- matrix(unlist(type1), nrow = length(type1), ncol = length(type1[[1]]), byrow = TRUE)
    dimnames(type1Mat) <- list(NULL, names(type1[[1]]))
    type1Mat <- type1Mat %>% as_tibble() %>% rowid_to_column("ID")
    type1Res <- survSimComb %>% left_join(type1Mat, by = "ID") %>% mutate(type1All = region_decline_mean) %>% dplyr::select(ID, type1All)

    # compile results
    resultMat <- matrix(unlist(result), nrow = length(result), ncol = length(result[[1]]), byrow = TRUE)
    dimnames(resultMat) <- list(NULL, names(result[[1]]))
    resultMat <- resultMat %>% as_tibble() %>% rowid_to_column("ID")
    finalRes <- survSimComb %>% left_join(resultMat, by = "ID") %>% mutate(monitBudget = budget, monitYears = monitRep - 1, survInter = monInter, monitStrat = strat, monitMethod = method, type2All = 1 - region_decline_mean,
                      type2GC = 1 - pop1_decline, type2KC = 1 - pop2_decline, type2MB = 1 - pop3_decline, type2NS = 1 - pop4_decline, type2N = 1 - pop5_decline, type2NW = 1 - pop6_decline,
                      type2SW = 1 - pop7_decline) %>%
                      left_join(type1Res, by = "ID") %>%
                      dplyr::select(monitBudget, monitYears, monitStrat, monitMethod, survIntens, survInter, type1All, type2All, type2GC, type2KC, type2MB, type2NS, type2N, type2NW, type2SW) %>%
                      as_tibble() %>% arrange(monitBudget, type2All)

    # save compiled results
    write.csv(x = finalRes, file = paste("output/results/ranking", "_decline", popSimComb[i, "decline"], "_envtlV",
                  popSimComb[i, "envtlV"], "_envtlCorr", popSimComb[i, "envtlCorr"], ".csv", sep = ""))
  }
}

# compile monitoring strategy details

# get surveys
Surveys <- apply(survSimComb, MARGIN = 1, FUN = function(x) {readRDS(file = paste("output/surveys/survey", "_budget", x["budget"], "_survIntens",
              as.numeric(x["survIntens"]), "_method", x["method"], "_monInter", as.numeric(x["monInter"]), "_strat", x["strat"], "_monitRep",
              as.numeric(x["monitRep"]), ".rds", sep = ""))})

# loop through monitoring strategies
for (i in 1:nrow(survSimComb)) {
  # get output
  Output <- Surveys[[i]] %>% mutate(HABITAT = Quality, LAND_USE = LUCategory, KOALA_DENS = KD_MEAN, NUM_SITES = SITES) %>%
                dplyr::select(FID, POP, HABITAT, LAND_USE, AREA, KOALA_DENS, NUM_SITES, SIZE, TYPE)

  # write output
  write.csv(x = Output, file = paste("output/surveys/monit_strat", "_monitBudget",
                  as.numeric(survSimComb[i, "budget"]), "_survIntens", as.numeric(survSimComb[i, "survIntens"]), "_monitMethod", survSimComb[i, "method"], "_survInter",
                  as.numeric(survSimComb[i, "monInter"]), "_monitStrat", survSimComb[i, "strat"], "_monitYears", as.numeric(survSimComb[i, "monitRep"]) - 1, ".csv", sep = ""))
}

# test adaptive strategy to abandon sites if koalas not seen for a number of surveys

# generate simulated survey data
# get survey information
# here we choose the survey strategies that usually come out as the top two based on Type 2 errors with a budget of 190 (see report) to evaluate the adaptive strategy
survChosen <- survSimComb[c(56, 65),]
Surveys <- apply(survChosen, MARGIN = 1, FUN = function(x) {readRDS(file = paste("output/surveys/survey", "_budget", x["budget"], "_survIntens",
              as.numeric(x["survIntens"]), "_method", x["method"], "_monInter", as.numeric(x["monInter"]), "_strat", x["strat"], "_monitRep",
              as.numeric(x["monitRep"]), ".rds", sep = ""))})
# loop through population combinations
for (i in 1:nrow(popSimComb)) {

  # get population simulations
  popSim <- readRDS(file = paste("output/popsims/popsim", "_decline", popSimComb[i, "decline"], "_envtlV",
        popSimComb[i, "envtlV"], "_envtlCorr", popSimComb[i, "envtlCorr"], ".rds", sep = ""))

  # generate data for the adaptive strategy
  # loop through adaptation years - here we use between 2 and 5
  for (yearAdapt in 2:5) {

    procTime <- system.time({

    cl <- makeCluster(detectCores() - 1)
    clusterExport(cl, list("simDataAdapt", "popSim", "f0", "f0se", "stripMissLow", "stripMissMean", "stripMissHigh", "yearAdapt"))
    clusterEvalQ(cl, library(tidyverse))
    clusterEvalQ(cl, library(extraDistr))
    clusterEvalQ(cl, library(unmarked))
    clusterEvalQ(cl, library(abind))
    dataTemp <- parLapply(cl = cl, X = Surveys, fun = function(x) {return(simDataAdapt(x, popSim, f0, f0se, stripMissLow, stripMissMean, stripMissHigh, yearAdapt))})
    stopCluster(cl)
    saveRDS(object = dataTemp, file = paste("output/adapt/data/simdata_adapt", "_decline", popSimComb[i, "decline"], "_envtlV",
          popSimComb[i, "envtlV"], "_envtlCorr", popSimComb[i, "envtlCorr"], "_yearAdapt", yearAdapt, ".rds", sep = ""))
    rm(dataTemp)
    gc()

    })

    print(procTime)

  }

  # generate data for the standard strategy
  procTime <- system.time({

  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, list("simData", "popSim", "f0", "f0se", "stripMissLow", "stripMissMean", "stripMissHigh"))
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(extraDistr))
  clusterEvalQ(cl, library(unmarked))
  dataTemp <- parLapply(cl = cl, X = Surveys, fun = function(x) {return(simData(x, popSim, f0, f0se, stripMissLow, stripMissMean, stripMissHigh))})
  stopCluster(cl)
  saveRDS(object = dataTemp, file = paste("output/adapt/data/simdata_stand", "_decline", popSimComb[i, "decline"], "_envtlV",
        popSimComb[i, "envtlV"], "_envtlCorr", popSimComb[i, "envtlCorr"], ".rds", sep = ""))
  rm(dataTemp)
  gc()

  })

  print(procTime)

  rm(popSim)
  gc()
}

# fit models to adaptive and standard strategies but only for the 0%, 10%, 25% and 50% declines
popChosen <- popSimComb[which((popSimComb[ ,"decline"] == 0) | (popSimComb[ ,"decline"] == 10) | (popSimComb[ ,"decline"] == 25) | (popSimComb[ ,"decline"] == 50)), ]

# loop through population combinations
for (i in 1:nrow(popChosen)) {

  # fit models for the adaptive strategy
  # loop through adaptation years - here we use between 2 and 10
  for (yearAdapt in 2:5) {

    procTime <- system.time({

    # get simulated data for adaptive strategy
    surveyData <- readRDS(file = paste("output/adapt/data/simdata_adapt", "_decline", popChosen[i, "decline"], "_envtlV",
          popChosen[i, "envtlV"], "_envtlCorr", popChosen[i, "envtlCorr"], "_yearAdapt", yearAdapt, ".rds", sep = ""))
    # PICK HOW MANY REPLICATES - edit out if using all replicates
    surveyData <- lapply(surveyData, FUN = function(x) {x[1:100]})
    expR <-  log((1 - (popChosen[i, "decline"] / 100)) ^ (1 / genYears))
    cl <- makeCluster(detectCores() - 1)
    clusterExport(cl, list("fitModels", "pops", "expR"))
    clusterEvalQ(cl, library(tidyverse))
    clusterEvalQ(cl, library(unmarked))
    clusterEvalQ(cl, library(binom))
    resultTemp <- parLapply(cl = cl, X = surveyData, fun = function(x) {return(fitModels(x, pops, expR))})
    stopCluster(cl)
    saveRDS(object = resultTemp, file = paste("output/adapt/results/result_adapt", "_decline", popChosen[i, "decline"], "_envtlV",
                  popChosen[i, "envtlV"], "_envtlCorr", popChosen[i, "envtlCorr"], "_yearAdapt", yearAdapt, ".rds", sep = ""))
    rm(surveyData)
    rm(resultTemp)
    gc()

    })

    print(procTime)
  }

  # fit models for the standard strategy
  procTime <- system.time({

  # get simulated data for adaptive strategy
  surveyData <- readRDS(file = paste("output/adapt/data/simdata_stand", "_decline", popChosen[i, "decline"], "_envtlV",
        popChosen[i, "envtlV"], "_envtlCorr", popChosen[i, "envtlCorr"], ".rds", sep = ""))
  # PICK HOW MANY REPLICATES - edit out if using all replicates
  surveyData <- lapply(surveyData, FUN = function(x) {x[1:100]})
  expR <-  log((1 - (popChosen[i, "decline"] / 100)) ^ (1 / genYears))
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, list("fitModels", "pops", "expR"))
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(unmarked))
  clusterEvalQ(cl, library(binom))
  resultTemp <- parLapply(cl = cl, X = surveyData, fun = function(x) {return(fitModels(x, pops, expR))})
  stopCluster(cl)
  saveRDS(object = resultTemp, file = paste("output/adapt/results/result_stand", "_decline", popChosen[i, "decline"], "_envtlV",
                popChosen[i, "envtlV"], "_envtlCorr", popChosen[i, "envtlCorr"], ".rds", sep = ""))
  rm(surveyData)
  rm(resultTemp)
  gc()

  })

  print(procTime)
}

# compile adaptive and standard monitoring strategy details for the 1st ranked and 2nd ranked strategy
compiledAdapt_Rank1 <- matrix(NA, nrow = 4, ncol = 6)
dimnames(compiledAdapt_Rank1) <- list(NULL, c("decline", "standard_error", "yearAdapt2_error", "yearAdapt3_error", "yearAdapt4_error",
                        "yearAdapt5_error"))
compiledAdapt_Rank1[ ,"decline"] <- c(0, 10, 20, 50)
compiledAdapt_Rank2 <- matrix(NA, nrow = 4, ncol = 6)
dimnames(compiledAdapt_Rank2) <- list(NULL, c("decline", "standard_error", "yearAdapt2_error", "yearAdapt3_error", "yearAdapt4_error",
                        "yearAdapt5_error"))
compiledAdapt_Rank2[ ,"decline"] <- c(0, 10, 20, 50)

# loop through declines
Declines <- c(0, 10, 25, 50)
for (i in 1:length(Declines)) {

  popChosenTemp <- popChosen[which(popChosen[, "decline"] == Declines[i]), ]

  #get the errors for the standard approach
  meanStand_Rank1 <- mean(apply(popChosenTemp, MARGIN = 1, FUN = function(x) {readRDS(file = paste("output/adapt/results/result_stand", "_decline",
              x["decline"], "_envtlV", x["envtlV"], "_envtlCorr", x["envtlCorr"], ".rds", sep = ""))[[1]]["region_decline_mean"]}))
  meanStand_Rank2 <- mean(apply(popChosenTemp, MARGIN = 1, FUN = function(x) {readRDS(file = paste("output/adapt/results/result_stand", "_decline",
                          x["decline"], "_envtlV", x["envtlV"], "_envtlCorr", x["envtlCorr"], ".rds", sep = ""))[[2]]["region_decline_mean"]}))
  if (i > 1) {
    meanStand_Rank1 <- 1 - meanStand_Rank1
    meanStand_Rank2 <- 1 - meanStand_Rank2
  }
  compiledAdapt_Rank1[i, "standard_error"] <- meanStand_Rank1
  compiledAdapt_Rank2[i, "standard_error"] <- meanStand_Rank2

  #get the errors for the adaptive approach
  for (j in 2:5) {
    meanAdapt_Rank1 <- mean(apply(popChosenTemp, MARGIN = 1, FUN = function(x) {readRDS(file = paste("output/adapt/results/result_adapt", "_decline",
                x["decline"], "_envtlV", x["envtlV"], "_envtlCorr", x["envtlCorr"], "_yearAdapt", j, ".rds", sep = ""))[[1]]["region_decline_mean"]}))
    meanAdapt_Rank2 <- mean(apply(popChosenTemp, MARGIN = 1, FUN = function(x) {readRDS(file = paste("output/adapt/results/result_adapt", "_decline",
                            x["decline"], "_envtlV", x["envtlV"], "_envtlCorr", x["envtlCorr"], "_yearAdapt", j, ".rds", sep = ""))[[2]]["region_decline_mean"]}))
    if (i > 1) {
      meanAdapt_Rank1 <- 1 - meanAdapt_Rank1
      meanAdapt_Rank2 <- 1 - meanAdapt_Rank2
    }
    compiledAdapt_Rank1[i, j + 1] <- meanAdapt_Rank1
    compiledAdapt_Rank2[i, j + 1] <- meanAdapt_Rank2
  }
}

# write outputs
write.csv(x = compiledAdapt_Rank1, file = "output/adapt/results/adapt_compiled_results_rank1.csv")
write.csv(x = compiledAdapt_Rank2, file = "output/adapt/results/adapt_compiled_results_rank2.csv")
