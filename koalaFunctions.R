# koalaFunctions 
# Code for simulating koala populations, surveys, and ...stats?  
# ... can maybe split into files in those subcategories later if this gets too big


require(extraDistr) # The "extraDistr" package is needed to use half-Normal distribution: surveyLines fn

#### Set some default values  #### 
doDefaultsHere <- 0 # This stops problems with hard-coding elsewhere in devpt
if (doDefaultsHere == 1){
    nRVs <- 1000.00 # no. RVs for this test
    dens <- 0.1 # density of koalas, per ha
    # strip sampling;
    areaStrip <- 70000   # matched to line transect params
    pObsStrip <- 0.84 # after Dique+al(2001)
    # line sampling;
    f0 <- 0.0285 # distance sampling parameter from Dique+al(2003); eff(W) = 1/f0
    f0se <- 0.0021 # S.E. for f0
    lenLine <- 1000 # 1km line transect length 
    # ... for f0 = 0.0285 this gives effW = 70m so effArea = 7 ha 
    
    # Pop trends
    dens0 <- 0.1  # density of koalas
    nYears <- 20 # years for simulation
    trend <- 0.9333 # 50% decrease in a decade
    demogV <- 1 # demographic stochasticity scaler
    envtlV <- 0.01 # environmental stochasticity scaler
}

surveyLines <- function(dens, lenLine, f0, reps=1){
  # Produce line transect observations - sightings and distances. 
  # Input dens (per HECTARE); lenLine (METRES); f0 (sightability param), reps (defaults to 1)
  
  # Output a list of nObsLine (#Observed per rep), and dObsLine (perpendicular distances for nObs>0)
  
halfNshape <- pi* f0^2 / 2 # 1 # precision (1/var) for halfNormal; tau in Rhodes et al. 2015 report
# /^\ The halfNormal mean is tau*sqrt(2/pi)

poisLam <- 2*lenLine*dens/(10000*f0)  
nObsLine <- rpois(n=reps, lambda=poisLam) # numbers observed
# /^\ Could possibly cope with DENS as a vector, but generating the list below would need edits

# Perpendicular distance of observations:
dObsLine <- rhnorm(n=sum(nObsLine), sigma=sqrt(1/halfNshape)) # distances observed: note draw not #Reps but the #observed

# Put distances in list for output (can count #obs from |dists| later)
startInd <- c(1, 1+cumsum(nObsLine)[1:reps-1]) # indices to split all the distances into their correct surveys
endInd <- cumsum(nObsLine)
lineLiszt <- vector("list", reps)
for (iRep in 1:reps){
  lineLiszt[iRep] <- list(dObsLine[ startInd[iRep] : endInd[iRep] ])
}

return(lineLiszt) # return the list of surveys' (vectors of) observed distances 
# fn formatDist.R changes these lists into Distance-ready tables 

} # end of fn surveyLines(dens, lenLine, f0)



surveyStrip <- function(dens, areaStrip, pObsStrip, reps=1){
# Produce strip transect observations #
popSizeStrip <- areaStrip*dens/10000

popSizeStrip <- round(popSizeStrip) #  for whole numbers of koalas 

nObStrip <- rbinom(n=reps, size=popSizeStrip, prob=pObsStrip)
return(nObStrip)
} # of fn surveyStrip(dens, areaStrip, pObsStrip, reps=1)

surveyAreas <- function(dens, areaCtArea, pObsArea, reps=1){
# simulate total-area count observations #
# same code as surveyStrip but option to enter different search-area and P(obs)
  popSizeArea <- areaCtArea*dens/10000
  
  popSizeArea <- round(popSizeArea) # whole koalas 
  
  nObsArea <- rbinom(n=reps, size=popSizeArea, prob=pObsArea)
  return(nObsArea)  
} # of fn surveyAreas(dens, areaCtArea, pObsArea, reps=1)

#### Produce population trajectories: 
projPops <- function(dens0, trend, nYears, demogV=0, envtlV=0)
  #  Input: dens0 in an array of #sites x #reps; and params trend, nYears, demographic & environmental variation
  # Output: all the population projections, as an array of size [#yrs+1 x #sites x #reps]
  #         to fit into a [3d-strata] x nLocs x y+1 x reps  e.g.  7x2x2x6x21x1000 = 3.5m elements
  {
  nLocs <-  max(1, dim(dens0)[1] ) # how many locations: max(1...) allows for single loc, scalar dens0 input
  nReps <-  max(1, dim(dens0)[2] ) # how many replicates 
  popArray <- array(NA, c(1+nYears, nLocs, nReps)) # set up results array
  popArray[1, , ] <- dens0  # initialise each loc/rep with its dens0
  # step through time now:
  for (tt in 1:nYears)
    {
    stdsNow <- sqrt(demogV /popArray[tt,,] + envtlV) # var = demV/N + envV
    stdsNow[is.infinite(stdsNow)] <- 0 # trap NaNs occurring when pop=0
    dPopRVs <- rnorm(n= nReps*nLocs, mean=trend, sd=stdsNow) # so R ~ N(R, var) like Lande 1993
    dPopsNow <- matrix( pmax(0,dPopRVs), nLocs, nReps) # reshape into matrix, limit new pops to >=0
    
    nextPops <- popArray[tt,,]*dPopsNow
    popArray[tt+1,,] <- nextPops
  } # end of looping for (tt in 1:nYears)
  
  popArray <- aperm(popArray, c(2,1,3)) # change shape from [yr x loc x rep] to [loc x yr x rep]
  dimnames(popArray) <- list("Location" = (1:nLocs),
                              "Year" = (1:(nYears+1)),
                              "Rep" = (1:nReps) )
  return(popArray)
} # end of fn projPops(dens0, trend, nYears, demogV=0, envtlV=0)



read.excel <- function(header=TRUE,...) {
  # read.excel() will paste stuff, with header(s), directly from Excel 
  read.table("clipboard",sep="\t",header=header,...)
}


# Spare code for ref, handling projPops output:
# # To rearrange into a yrs x reps (ignore location?): maybe use as.tibble (incl with as.tbl_cube):
# # https://community.rstudio.com/t/fastest-way-to-convert-3d-array-to-matrix-or-data-frame/38398/2
# # dimnames(popArray) <- list("Yrs" = sprintf("Year%d", 1:(nYears+1)),
# #                                    "Reps" = sprintf("Rep%d", 1:nReps),
# #                                    "Locs" = sprintf("Loc%d", 1:nLocs))
# 
# 
# popTibcub <- as.tbl_cube(popArray)
# popTibble <- as_tibble(popTibcub)
# # then Something (not quite) like: 
# ggplot(popTibble, aes(x=Year, y=popArray), group_by(Location)) + geom_line()
# 
# # lines(x=matrix(rep(0:nYears, nReps), nrow=nReps), y=popArray[,,1])
# 
# library(ggthemes) # fancy themes - Tufte, 385 etc.
# library(plotly) # interactive graphics
