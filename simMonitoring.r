# load required libraries
library(tidyverse)
library(unmarked)
library(MASS)
library(parallel)

# location of functions
#source("functions.r")

# define input parameters
genYears <- 20 # number of years for a generation
demogV <- 96.21 # demographic stochasticity variance - this is variance in individual annual growth rate among individuals - taken from Rhodes et al. (2011) based on the radio tacking data only
                # model. Here the standard deviation of the growth rate parameter estimate was 0.04458435 from 220 individuals, so we calculated demogV as (0.04458435 * 220) ^ 2
envtlV <- c(0, 0.01, 0.02) # environmental stochasticity variance
envtlCorr <- c(0, 0.5, 1) # correlation in environmental stochasticity among populations
nSimsPop <- 10000 # number of replicates for the population simulations
nYears <- 50 # simulation time horizon
decline <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50) # percentage decline over three generations
f0 <- 0.0285 # distance detection parameter from Dique et al(2003), eff(W) = 1/f0
f0se <- 0.0021 # standard error of f0 from Dique et al(2003)
tranDayLine <- 4600 # line transect length per day in m
tranDayArea <- 7.5 # strip transect / area search area per day in ha
siteSize <- 100 # site size in ha
monitRep <- c(5, 10, 15, 20) # monitoring reporting years
budget <- c(95, 190, 380) # field days per year
survIntens <- c(30, 50, 70) # proportion of the site surveyed
method <- c("line", "area") # "line" is all line and "area" is all area searches
monInter <- c(1, 2) # monitoring interval in years
strat <- c("equal", "area", "density", "invdensity")

# load strata and initial population densities - set -9999 to NA and remove strata with < 100 ha of habitat
# then group by genetic population and nest
# note here that we replace the density data for "core - urban footprint - north west" stratum with the density data for "core - urban footprint - north west" stratum
# since this density data is missing for "core - urban footprint - north west". we also use data from Biolink 2019 for North Stradbroke Island (NSI) and assume core and
# non-core habitat densities are the same for NSI. this is the file "start_dens_modified.csv"
pops <- read.csv("input/start_dens_modified.csv") %>%
  mutate(KD_MEAN=replace(KD_MEAN, KD_MEAN==-9999, NA), KD_SD=replace(KD_SD, KD_SD==-9999, NA),
    KDLN_MEAN=replace(KDLN_MEAN, KDLN_MEAN==-9999, NA), KDLN_SD=replace(KDLN_SD, KDLN_SD==-9999, NA)) %>%
  filter(AREA >= 100) %>% mutate(ABUND = round(KD_MEAN * AREA)) %>% as_tibble() %>% rowid_to_column("ID") %>% group_by(POP) %>% nest()
# calculate the population size for each genetic populations
pops <- pops %>% mutate(N = map_dbl(data, ~sum(.$ABUND)))

# simulate population dynamics - comment out if already run
popSimComb <- expand.grid(decline = decline, envtlV = envtlV, envtlCorr = envtlCorr)
#cl <- makeCluster(detectCores() - 1)
#clusterExport(cl, list("pops", "projPops", "nSimsPop", "nYears", "genYears", "demogV" ))
#clusterEvalQ(cl, library(tidyverse))
#clusterEvalQ(cl, library(MASS))
#parApply(cl = cl, X = popSimComb, MARGIN = 1, FUN = function(x) {saveRDS(object = projPops(pops = pops$N, nYears = nYears, decline = x["decline"], genYears = genYears,
#                demogV = demogV, envtlV = x["envtlV"], envtlCorr = x["envtlCorr"], nSims = nSimsPop), file = paste("output/popsims/popsim", "_decline",
#                x["decline"], "_envtlV", x["envtlV"], "_envtlCorr", x["envtlCorr"], ".rds", sep = ""))})
#stopCluster(cl)

# generate surveys - comment out if already run
# get combinations of parameters
survSimComb <- expand.grid(budget = budget, survIntens = survIntens, method = method, monInter = monInter, strat = strat, monitRep = monitRep)
# collate fixed survey parameters together
fixedSurvParams <- c(f0, f0se, tranDayLine, tranDayArea, siteSize)
names(fixedSurvParams) <- c("f0", "f0se", "tranDayLine", "tranDayArea", "siteSize")
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, list("getSurveys", "smart.round", "pops", "survSimComb", "fixedSurvParams" ))
clusterEvalQ(cl, library(tidyverse))
parApply(cl = cl, X = survSimComb, MARGIN = 1, FUN = function(x) {saveRDS(object = getSurveys(pops = pops, survParams = c(x, fixedSurvParams)), file = paste("output/surveys/survey", "_budget",
                x["budget"], "_survIntens", as.numeric(x["survIntens"]), "_method", x["method"], "_monInter", as.numeric(x["monInter"]), "_strat", x["strat"], "_monitRep", as.numeric(x["monitRep"]), ".rds", sep = ""))})
stopCluster(cl)
