# demoKoala
# Demo of sampling strategies: Rotate through 3 zones vs. Visit all zones annually
# Code runs OK if two components are tweaked/skipped;
# 1. "NA" for South population's urban core habitat initial density (tweak using Lines 91-92 avDens0[8] etc. ) 
# 2. Distance observations should be combined, too few prevent model fitting by Distance package so code aborts
# ... if [1] is tweaked (assign false initial densities for 8th stratum) 
# ... and [2] skipped (comment out lines c.200); get 6170 observations overall

set.seed(42)

require(Distance)
# Get input/starting variables
trend <- 0.5^(1/10) # = 0.933; annual trend, e.g. 50% decrease in a decade
nYears <- 9
f0 <- 0.0285    # distance detection parameter
demogV <- 0     # demographic stochasticity scaler - difficult with densities...
envtlV <- 0.01  # environmental stochasticity scaler

arbitrarySiteArea <- 10*1e6 # m^2 so (A x 1e6) = A sq.km

# *****************************************************************************
#### Establish monitoring schedules ####
# *****************************************************************************
# Year - PopLocn - ifCore - ifUrban - Site# - Effort - (Method)
# [ 20  x   3       x 2         x 2    (1-3/1-9)   (varies) (Line) ]
# EXAMPLE:
# [0,   "North",   "Core",  "Bush", 1, 5000[m], ("Line")]
# [0,   "North",   "Core",  "Bush", 2, 5000[m], ("Line")]
# [ .... etc.]
# [1,   "North",   "Core",  "Bush", 1, 5000[m], ("Line")]
# [1,   "North",   "Core",  "Bush", 2, 5000[m], ("Line")]
# [ .... etc.]
# IGNORING FOR NOW (development): 
#  within-year timesteps (need to alter pop model for that)
#  strip and area sampling (although they're coded up - see koalaFunctions.R)

print("Enter T/F here to swap plans:")
FewSitesLots <- FALSE # Either [0] do 9 sites in 1 popZone a year, revisiting triennially; or [1] 3 sites in 3 popZones annually

sampYrs  <- 1:nYears
PopZones <- c("North","South","West")
ifCore <- c("Core", "Noncore")
ifBush <- c("Bush", "Urban")
survLocs <- "set in IF statement below"
Efforts <- 5000 # line transect length in metres
Methods <- "Line"
nReps <- 10

if (FewSitesLots){ # if visiting a few sites in each Zone each year
    survLocs <- 1:3
    survPlan <- expand.grid(sampYr=sampYrs, PopZone=PopZones, ifCore=ifCore, ifBush=ifBush, 
                        surLocNo = survLocs, Effort=Efforts, trueDens = NA)#, Method = Methods
                # ^^ could also convert to official tibble: as.tbl(survPlan)
} else { # if rotating through pop-zones
    survLocs <- 1:9
    survPlan <- expand.grid(sampYr=sampYrs, PopZone="FillInBelow", ifCore=ifCore, ifBush=ifBush, 
                            surLocNo = survLocs, Effort=Efforts, trueDens = NA)#, Method = Methods

    levels(survPlan$PopZone) <- PopZones # or? c(levels(survPlan$PopZone), PopZones)     
    # Now loop through years assigning zones:
    for (yy in sampYrs){
        zoneIndex <- length(PopZones) - (-yy %% length(PopZones)) # e.g. will give 1 2 3 1 2 3....  as loop progresses

        survPlan$PopZone[survPlan$sampYr==yy  ] <- PopZones[zoneIndex] # Assign this year's zone 
    }
}
nSurvTypes <- dim(survPlan)[1]

#   Survey schedule done... but can replace Effort with +/-effortRV maybe, or Line with Strip for non-core etc.


# *****************************************************************************
#### POPULATIONS ####
# ***************************************************************************** 
# Run population model to get population densities w.r.t. monitoring schedule 
# Here have 12 strata x (nLocs=3/9) sites over 9 years; and RR replicates
# Assuming max once-a-year monitoring.

nLocs <- max(survLocs)

# *****************************************************************************
#### Initial Densities #### 
# *****************************************************************************
# Generate initial pop density for each location and rep:
for (tempSetupWillReadFromFileIdeally in 1:1){
    
    avDens0 <- c(0.038769138, 0.047350915, 0.031528022, 0.03070217, # North (Moreton Bay)
             0.027137499, 0.041346916, 0.158925215, NA,        # South (Gold Coast)
            0.020715351, 0.02314302, 0.050048845, 0.033880637) # West (SW)
    sdDens0 <- c(0.097636339, 0.03772673, 0.050980236, 0.011581464,
                0.067281493, 0.037211536, 0.49413686, NA,
                0.019543029, 0.015604415, 0.101918524, 0.010924325)
# avDens0[8] <- 0.05 # overcome NA hassle
# sdDens0[8] <- 0.10 # overcome NA hassle

            
    # Need dens0 for strata x locations 
    # projn fn will produce yrs [x location] x reps
    dens0stats <- expand.grid(ifBush=ifBush, ifCore=ifCore, PopZone=PopZones, 
                             means = NA, stdvns = NA)   # Make factorial grid of entries (PopZones last arg => adjacent)
    dens0stats <- dens0stats[ , c(3,2,1,4,5)]             # Swap order to put PopZones first etc.
    
    dens0stats$locParam <- log(avDens0^2 / sqrt(sdDens0^2 + avDens0^2)) # R-friendly parameters for logN distn
    dens0stats$shaParam <- sqrt(log(1 + (sdDens0^2 / avDens0^2)))  # its "location" parameter & shape parameter
    # see https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/  (and Burgman+al'93)
    
} # end of "temp-Setup-Will-Read-From-File-Ideally"


# *****************************************************************************
#### Population Projection #### 
# *****************************************************************************
# Calls my function "projPops"

# ptm <- proc.time()  # TIMING : Start the clock

densArray <- array(NA, c(length(PopZones), length(ifCore), length(ifBush), nLocs, nYears+1, nReps))
# set up array of size [pops x hab x urbush x locs x yrs x Reps] e.g. [7 x 2 x 2 x 6 x 21 x 1000] 
# ... and now fill it in with simulated popualtion densities:
for (iZone in 1:length(PopZones)){
    zoneNow <- PopZones[iZone]
    for (iCore in 1:length(ifCore)){
        coreNow <- ifCore[iCore]
        for (iBush in 1:length(ifBush)){
            bushNow <- ifBush[iBush]
            stratumNow <- c(zoneNow, coreNow, bushNow) # just note stratum for quick ref if needed
            tempSlice <- dens0stats %>% filter(PopZone==zoneNow, ifCore==coreNow, ifBush==bushNow)
            # Get starting densities by drawing nSites x nReps from relevant distribution
            stratLNm <- tempSlice$locParam
            stratLNs <- tempSlice$shaParam
            D0mat <- rlnorm(n= nLocs*nReps, meanlog=stratLNm, sdlog=stratLNs)  #  draw initial densities; to make matrix:
            D0mat <- matrix(D0mat, nLocs, nReps) # reshape into [nLocs x nReps] matrix of starting densities
            # did use max0 - not sure why? D0mat <- matrix( pmax(0,D0mat), nLocs, nReps) # reshape to [nLocs x nReps].. 
            
            popProj <- projPops(D0mat, trend, nYears, demogV, envtlV) # project D0 matrix that by years and replicates 
            densArray[iZone, iCore, iBush, , , ] <- popProj # put popProj [#Locs x #Yrs+1 x #Reps] array into its stratum
        } # end (non)Bush loop
    } # end (non)core loop
} # end PopZone loop
handyAv <- rowMeans(popProj, dims=2) # handy to see mean trajectory for each location in last stratum calc'd above.
handySD <- apply(popProj, c(1,2), sd) # (just noting alternative fn method to get popProj stats, like SD here)
rm(popProj) # can be big, e.g. 10MB for 6 sites x 21 yrs x 10000 reps; then densArray is [~that x 28] 


# Now have 
# survPlan ~ "sampYr"   "PopZone"  "ifCore"   "ifBush"   "surLocNo" "Effort"   "trueDens"(=NA still)
# and densArray ~  zones x (non)Core x (non)Urban  x Locs x  Yrs+1 x Reps
#
# So can do survey inputting true density ... and recording estimated density
#-OR- [see end / koalaMain for alternative approach]


# *****************************************************************************
#### Go Monitoring #### 
# *****************************************************************************

# Step through monitoring schedule to generate survey results
# calls my function surveyLines() to generate surveys
# calls my function formatDist() to reshape those results for the Distance package
# and then calls fn Distance::ds() to do distance analysis

totalSeen <- 0 # running total of observations, to check if it's working at all!


rowNA <- rep( NA, nReps*nrow(survPlan) ) # empty vector
estDensTab <- data.frame("survYr" = rowNA, "PopZone" = rowNA, "ifCore" = rowNA, "ifBush" = rowNA, "locNo" = rowNA, 
                 "Rep" = rowNA, "ActDens" = rowNA, "EstDens" = rowNA) # empty holder for all survey results
 
for (ss in 1:nrow(survPlan)){
    thisSurvey <- survPlan[ss,]
    
    indZone <- match(thisSurvey$PopZone, PopZones) # index into c("North","South","West")
    indCore <- match(thisSurvey$ifCore, ifCore)
    indBsUr <- match(thisSurvey$ifBush, ifBush)
    indYr <- thisSurvey$sampYr +1 # add one because first index in densArray is for year=0

    # densArray structure ~  [zones x (non)Core x (non)Urban  x Locs x  Yrs+1 x Reps], so...
    repDens <- densArray [ indZone, indCore, indBsUr , thisSurvey$surLocNo, indYr ,   ] # nReps # of true densities
    lenLine = thisSurvey$Effort # length of line transect for this survey
    
    resIndex <- (1:nReps) + (ss-1)*nReps # indices to record this survey-types results
    estDensTab[resIndex, "survYr"] <- thisSurvey$sampYr # set up results with survey type info
    estDensTab[resIndex, "PopZone"] <- thisSurvey$PopZone 
    estDensTab[resIndex, "ifCore"] <- thisSurvey$ifCore 
    estDensTab[resIndex, "ifBush"] <- thisSurvey$ifBush
    estDensTab[resIndex, "locNo"] <- thisSurvey$surLocNo 
    
    for (rr in 1:nReps){  
        survResult <- surveyLines(repDens[rr], lenLine, f0, reps=1) # do a line transect ...:
        # Outputs a list of perpendicular distances in m^2 (for nObs>0; otherwise NA)
        
        totalSeen <-totalSeen + length(unlist(survResult)) # running total of no.observations (for checking) 
        
        survRes4dist <- formatDist(survResult, locArea = arbitrarySiteArea, Effort = lenLine) 
        # ^^ re-format the results for DISTANCE package

        
        #
        # NB .... This bit may be done prematurely; 
        # (may need to wait and combine all results (for year?) for distance model to fit)
        #
        distOut <- ds(data=survRes4dist, key="hn", adjustment=NULL) # DISTANCE::ds(); fit with unadjusted half-Normal
        summOut <- summary(distOut)

        densum <- summOut$dht$individuals$D # density output, e.g. (minke whale example from vignette):
        #   Label   Estimate          se        cv         lcl        ucl       df
        # 3 Total 0.02403400 0.007179465 0.2987212 0.012838347 0.04499280 14.00459
        #= (Stratum)  ^den^      s.e.        c.v.    Lower CI     Upper CI   d.f.

        # Fill in results:
        resIndex <- rr + (ss-1)*nReps
        estDensTab[resIndex, "Rep"] <- rr # record rep no.
        estDensTab[resIndex, "ActDens"] <- repDens[rr] # record true density
        estDensTab[resIndex, "EstDens"] <- densum$Estimate # record est. density
        
    } # end loop: for (rr in 1:length(repDens))
} # end loop: for (ss in 1:nrow(survPlan))

# NB checking units:
 print("dens is /ha above, surveyLines() corrects to m2 so Distance::ds() uses that; ensure site area in m^2 too so")
# dist_units <- print("all done in metres here; see vignette's wren_lt_units if nec.")  

 
# Calculate trend with stats model
 # ...
# Detect change or not - record this
 # ...
# Find power
 # ...

  
#### Old Spare code #### 
 
# # Earlier ln-density set ups: (generated with   exp( rnorm(n= nLocs*nReps, mean=stratAv, sd=stratSD))  ) 
 # avsLnD <- c(-5.3, -5.4, -5.5,-5.6,-5.4,-5.5,-5.6,-5.7,-3.6,-3.7,-3.8,-3.9) # logD is ~ N(avsLnD,sdsLnD)
 # sdsLnD <- c(2.1,2.0,1.9,1.8,1.3,1.2,1.1,1.0,1.0,0.9,0.8,0.7) # made these up in koalaCalcs/simDens
 
 # # Or generate random log-densities:
 # avsLnD <- runif(4*length(PopZones), min= -0.7, max = -0.2) # [-0.7,-0.2] -> dens=[0.5,0.8]
 # sdsLnD <- runif(length(avsLnD), min=0.5, max = 2.5)
 # lnD0stats <- expand.grid(ifBush=ifBush, ifCore=ifCore, PopZone=PopZones, 
 #                          means = NA, stdvns = NA)   # Entry order here reflected vector above from koalaCalcs/simDens
 # lnD0stats <- lnD0stats[ , c(3,2,1,4,5)]             # So swap order of columns to match demo stats list ;-/
 # lnD0stats$means <- avsLnD
 # lnD0stats$stdvns <- sdsLnD
 # #
 # print("I just made up arbitrary lnDensities - maybe reset with strata-fied ones for more realism")
 
 
 ## This code visits each site/replicate once a year:
 ## Old, started looping through strata:
 # for (iPopn in 1:length(PopZones)){     thisPop <- PopZones[iPopn]
 #     for (iCore in 1:length(ifCore)){        thisCnC <- ifCore[iCore]
 #         for (iBush in 1:length(ifBush)){            thisBvU <- ifBush[iBush]
 #             thisStrat <- c(thisPop, thisCnC, thisBvU) # just for quick ref if needed
 #             stratumPlan <- survPlan %>% filter(PopZone==thisPop, ifCore==thisCnC, ifBush == thisBvU)
 #         } # end: for (iBush in 1:length(ifBush))
 #     } # end: for (iCore in 1:length(ifCore))
 # } # end: for (iPopn in 1:length(PopZones))
 
 