# koalaMain
set.seed(42)
# Shell script to call other functions.  In devpt.  "demoKoala" is more developed, for 2 selected monitoring strategies and pared-down .

# IGNORING FOR NOW (development): 
#  within-year timesteps (need to alter pop model for that)
#  strip and area sampling (although they're coded up - see koalaFunctions.R)

# Get input/starting variables
trend <- 0.5^(1/10) # = 0.933; annual trend
nYears <- 5

demogV <- 1 # demographic stochasticity scaler
envtlV <- 0.01 # environmental stochasticity scaler


#### Establish *one* monitoring schedule (or read/write a file for it): ####
# Year - subYr - PopLocn - ifCore - ifUrban - Site# - Effort - (Method)
# [ 20  x  2(?)  x   7       x 2         x 2    (?)   (varies) (Line) ]
# EXAMPLE:
# [0,   (7),  "Noosa",   "Core",  "Bush", 1, 500[m], ("Line")]
# [0,   (7),  "Noosa",   "Core",  "Bush", 2, 500[m], ("Line")]
# [1,   (1),  "Noosa",   "Core",  "Bush", 1, 500[m], ("Line")]
# [1,   (1),  "Noosa",   "Core",  "Bush", 2, 500[m], ("Line")]
# [ .... etc.]; and maybe
# [1,   (1),  "Noosa",   "Core",  "Urban", 2, 5[m2], ("tArea")]
# [ .... etc.]
sampYrs  <- 1:nYears
PopZones <- c("Noosa","GoldCoast")
ifCore <- c("Core", "Noncore")
ifBush <- c("Bush", "Urban")
survLocs <- 1:2
Efforts <- 50
Methods <- "Line"
survPlan <- expand.grid(sampYr=sampYrs, PopZone=PopZones, ifCore=ifCore, ifBush=ifBush, 
                        surLocNo = survLocs, Effort=Efforts, trueDens = NA)#, Method = Methods
# could also convert to official tibble: as.tbl(survPlan)
nSurvTypes <- dim(survPlan)[1]


#... and now can replace Effort with +/-effortRV maybe, or Line with Strip for non-core etc.

#### POPULATIONS - DONE ####
# Run population model to get population densities w.r.t. monitoring schedule ?

# Will have 28 strata x 6? sites over 20? years; and 1000 replicates?
# Where 6 is the *max* number of sites visited in *any* of the monitoring strategies
# And assuming max once-a-year monitoring.
# Even if e.g. Olympic-year surveys, have to calculate each year's pop to calculate demog stoch(?)
# Storage array is 28*6*21*1000 = 3,528,000

nLocs <- 6 # sites per stratum to visit (for Dev; will change with survey design)
nReps <- 100 # replicates (for Dev)

# Generate initial pop density for each location and rep:
for (tempSetupWillReadFromFileReally in 1:1){
    PopZones <- c("Gold","Noosa","FiveMore")
    ifCore <- c("Core", "Noncore")
    ifBush <- c("Bush", "Urban")
    avsLnD <- c(-5.3, -5.4, -5.5,-5.6,-5.4,-5.5,-5.6,-5.7,-3.6,-3.7,-3.8,-3.9) # logD is ~ N(avsLnD,sdsLnD)
    sdsLnD <- c(2.1,2.0,1.9,1.8,1.3,1.2,1.1,1.0,1.0,0.9,0.8,0.7) # made these up in koalaCalcs/simDens
    
    PopZones <- c("GoldCoast","KoalaCoast", 'MoretonBay', "Noosa","SouthWest","NorthWest","Straddie")
    
    avsLnD <- runif(4*length(PopZones), min= -6, max = -3.5)
    sdsLnD <- runif(length(avsLnD), min=0.5, max = 2.5)
    #
    #
    print("I just made up arbitrary lnDensities")
    #
    # Need dens0 for strata x locations 
    # projn fn will produce yrs [x location] x reps
    lnD0stats <- expand.grid(ifBush=ifBush, ifCore=ifCore, PopZone=PopZones, 
                             means = NA, stdvns = NA)
    lnD0stats <- lnD0stats[ , c(3,2,1,4,5)] # Swap order of columns to match demo stats list ;-/
    lnD0stats$means <- avsLnD
    lnD0stats$stdvns <- sdsLnD
} # end of "temp-Setup-Will-Read-From-File-Really"


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
            stratumNow <- c(zoneNow, coreNow, bushNow) # just for quick ref if needed
            tempSlice <- lnD0stats %>% filter(PopZone==zoneNow, ifCore==coreNow, ifBush==bushNow)
            # Get starting densities by drawing nSites x nReps from relevant distribution
            stratAv <- tempSlice$means
            stratSD <- tempSlice$stdvns
            D0mat <- exp(rnorm(n= nLocs*nReps, mean=stratAv, sd=stratSD))  #  draw starting densities; will become matrix ...
            D0mat <- matrix( pmax(0,D0mat), nLocs, nReps) # reshape into [nLocs x nReps] matrix of starting densities
            
            popProj <- projPops(D0mat, trend, nYears, demogV, envtlV) # project that by years and replicates 
            densArray[iZone, iCore, iBush, , , ] <- popProj # put popProj [#Locs x #Yrs+1 x #Reps] array into its stratum
        } # end (non)Bush loop
    } # end (non)core loop
} # end PopZone loop
handyAv <- rowMeans(popProj, dims=2) # see the mean trajectory for each location in last stratum calc'd above.  
handySD <- apply(popProj, c(1,2), sd) # Alternative setup to get stats 
rm(popProj) # can be big, e.g. 10MB for 6 sites x 21 yrs x 10000 reps; then densArray is [~that x 28] 

# proc.time() - ptm # TIMING : Stop+output the clock



#-OR- do within each surveyPlan?
# Run population model to get population densities w.r.t. monitoring schedule:
# Currently [3.8.2020] per stratum as it'll cope with multiple locations within a stratum
# Loop through strata:
# for (iPopn in 1:length(PopZones)){
#   thisPop <- PopZones[iPopn]
#   for (iCore in 1:length(ifCore)){
#     thisCnC <- ifCore[iCore]
#     for (iBush in 1:length(ifBush)){
#       thisBvU <- ifBush[iBush]``
#       stratumPlan <- survPlan %>% filter(PopZone==thisPop, ifCore==thisCnC, ifBush == thisBvU)
#       nLocsStrat <- max(stratumPlan$surLocNo) # # survey locations for which we need to ...
#       print("Something like: 
#             densStrat <- popmodel(dens0, trend, nYears, nLocsStrat, reps=1)")  # generate density trajecctories of 1:nYears
#       print("But may be better to do separately during set-up:  
#             (if doing 1000 reps, 28 strata, 6 sites, 20 years 
#             => 3,360,000 RV draws)")
#     } # end: for (iBush in 1:length(ifBush))
#   } # end: for (iCore in 1:length(ifCore))
# } # end: for (iPopn in 1:length(PopZones))


# Step through monitoring schedule to generate survey results
print("Ready to do this next... survey code already done in koalaFunctions")
# survResult <- surveyLines(dens, lenLine, f0, reps=1)

# This code visits each site/replicate once a year:
# Loop through strata:
for (iPopn in 1:length(PopZones)){
    thisPop <- PopZones[iPopn]
    for (iCore in 1:length(ifCore)){
        thisCnC <- ifCore[iCore]
        for (iBush in 1:length(ifBush)){
            thisBvU <- ifBush[iBush]
            thisStrat <- c(thisPop, thisCnC, thisBvU) # just for quick ref if needed
            stratumPlan <- survPlan %>% filter(PopZone==thisPop, ifCore==thisCnC, ifBush == thisBvU)
            nLocsStrat <- max(stratumPlan$surLocNo) # # survey locations for which we need to ...
            print("Something like:
            densStrat <- popmodel(dens0, trend, nYears, nLocsStrat, reps=1)")  # generate density trajecctories of 1:nYears
            print("But may be better to do separately during set-up:
            (if doing 1000 reps, 28 strata, 6 sites, 20 years
            => 3,360,000 RV draws)")
        } # end: for (iBush in 1:length(ifBush))
    } # end: for (iCore in 1:length(ifCore))
} # end: for (iPopn in 1:length(PopZones))


# Feed survey results into stats model



# Detect change or not - record this

# Find power.