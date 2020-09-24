formatDist <- function(transLiszt, locArea = -1, Effort = -1)
# formatDist(transLiszt, locArea = -1, Effort = -1)
# wrangle simulated koala line-survey data into format for use in DISTANCE package
# "transLizst" is output from surveyLines(), lists of sighting distances    
# "locArea" should be in m^2, Effort in metres.
  { 
# Help and source at end of file... 
# e.g. to test:
    # stestX <- surveyLines(1,5000,0.028,10)
    # transLiszt <- stestX 
    # locArea <- 42
    # print("hardcoded arbitrary lists and areas in formatDist")
    

    # Get properties of input lists:
    nLists <- length(transLiszt) # number of surveys
    nAllObs <- length(unlist(transLiszt)) # total # obsvns
    
    # Set up output dataframe "dataForDISTANCE"
    Region.Label <- rep("irrelvt?", nAllObs) # vector of stratum name
    Area        <- rep(locArea,  nAllObs) # vector of this site's area
    Sample.Label <- rep(-99,  nAllObs) # vector of transect numbers
    Effort      <- rep(lenLine,  nAllObs) # vector of effort (line length in metres)
    object      <- rep(-99,  nAllObs) # vector of sighting record number
    distance    <- rep(-99,  nAllObs) # vector of perpendicular distances, to fill in below
    Study.Area  <- rep("irrelvt?",  nAllObs) # vector of site names
    
    dataForDISTANCE <- data.frame(Region.Label, Area, Sample.Label, Effort, object, distance, Study.Area)
    
    # loop through obsvns and fill them in to dataForDISTANCE columns: Sample.Label, object and distance
    latestRow <- 0 # to count sequential rows as dataForDISTANCE is populated
    for(tt in 1:nLists){ 
        thisList <- unlist(transLiszt[tt])
        nObsHere <- length(thisList)
        rowInds <- latestRow+(1:nObsHere)
        
        dataForDISTANCE[rowInds, "Sample.Label"] <- tt # vector of this transect number
        dataForDISTANCE[rowInds, "object"]      <- 1:nObsHere # vector of sighting record numbers
        dataForDISTANCE[rowInds, "distance"]    <- thisList # vector of perpendicular distances, to fill in below
        
        latestRow <- max(rowInds) # set up for next list
    } # end of for(tt in 1:nLists) looping through transects  
    return(dataForDISTANCE)
    
}  # End of fn: formatDist(transLiszt)

# quoting the Distance package paper, Miller et al. (2019):
# To include additional covariates into the detection function a data.frame is required. 
#  Each row of the data.frame contains the data on one observation.
# The data.frame must contain a column named 
# DISTANCE (containing the observed distances) 
#  and additional named columns for any covariates that may affect detectability (for example observer or seastate). 
# [Reserved/unnecessary column names to avoid here: SIZE, OBJECT and DETECTED].
# 
# To estimate density or to estimate abundance beyond the sampled area, 
#  additional columns should be included in the data.frame specifying: 
# SAMPLE.LABEL, the ID of the transect;                   <- my Location-within-stratum?
# EFFORT, transect effort (for lines, their length);      <- lineLength
# REGION.LABEL, the stratum containing the transect;      <- stratum
# AREA, the area of the strata.                           <- Hmmmmm.... check jargon consistency
# 
# Transects which were surveyed but have no observations must be included 
#    in the data set with NA for distance and any other covariates. 
#... 
# If distances from a line  transect survey are recorded in meters, 
# the Effort columns should contain line lengths also in meters 
# and the Area column gives the stratum areas in square meters. 
# This would lead to density estimates of animals per square meter.

# > head(wren_lt); wren_lt[15*(2:4),]
# #      Region.Label    Area  Sample.Label Effort object distance Study.Area
# #   1     Montrave     33.2            1  0.416      5       15 Montrave 4
# #   2     Montrave     33.2            1  0.416      6       80 Montrave 4
# #___3_____Montrave     33.2            1  0.416      7       35 Montrave 4___
# # 30      Montrave     33.2            3  0.802     75        0 Montrave 4
# # 45      Montrave     33.2            5  0.700    115       35 Montrave 4
# # 60      Montrave     33.2            6  0.802    150       30 Montrave 4
