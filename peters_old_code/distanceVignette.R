# distanceVignette
# code entered from https://examples.distancesampling.org/Distance-lines/lines-distill.html
#  and from Miller et al. 2019 too
require("Distance")

data("wren_lt") #  LT stands for line transect ...
# > str(wren_lt)
# 'data.frame':	156 obs. of  7 variables

# > str(wren_lt_units)
# 'data.frame':	3 obs. of  3 variables: Variable  : chr  "Area" "Effort" "distance"
# > wren_lt_units
# Variable     Units Conversion
# 1     Area   hectare      10000
# 2   Effort kilometer       1000
# 3 distance     meter          1

head(wren_lt); wren_lt[15*(2:8),]; tail(wren_lt)
# Region.Label    Area  Sample.Label Effort object distance Study.Area
#   1     Montrave 33.2            1  0.416      5       15 Montrave 4
#   2     Montrave 33.2            1  0.416      6       80 Montrave 4
#___3_____Montrave 33.2            1  0.416      7       35 Montrave 4___
# 30      Montrave 33.2            3  0.802     75        0 Montrave 4
# 45      Montrave 33.2            5  0.700    115       35 Montrave 4
# 60      Montrave 33.2            6  0.802    150       30 Montrave 4
# 75      Montrave 33.2            7  0.786    178       80 Montrave 4
# 90      Montrave 33.2            8  0.810    205       65 Montrave 4
# 105     Montrave 33.2           12  0.094    237       90 Montrave 4
#_120_____Montrave 33.2           14  0.542    267       35 Montrave 4___
# 153     Montrave 33.2           18   0.40    338       50 Montrave 4
# 154     Montrave 33.2           19   0.04    340       40 Montrave 4
# 155     Montrave 33.2           19   0.04    341       80 Montrave 4
# 156     Montrave 33.2           19   0.04    343       25 Montrave 4

hist(wren_lt$distance, xlab="Distance (m)", main="Winter wren line transects")
# Seems bad obsvn on the line... evidence of evasive behaviour

# NB
conversion.factor <- convert_units(distance_units = "metre", effort_units = "kilometre", area_units = "hectare")
# NB


wren.hn <- ds(data=wren_lt, key="hn", adjustment=NULL, # HN=half-Normal; alternatively HR (=hazard rate), or unif[orm] 
              convert.units=conversion.factor)
# Fitting half-normal key function
# Key only model: not constraining for monotonicity.
# AIC= 1418.188
#
# > wren.hn
# Distance sampling analysis object
# Detection function:   Half-normal key function 
# Estimated abundance in covered region: 227.7249

# See webpage for more (not very helpful) exposition

# "Visually inspect the fitted detection function with the plot() function, 
# specifying the cutpoints histogram with argument breaks":
  cutpoints <- c(0,5,10,15,20,30,40,50,65,80,100)
plot(wren.hn, breaks=cutpoints, main="Half normal model, wren line transects")

# Note to self: 10 runs of each took 3.85s (halfN), 32.9s (Haz), and 32.2s (Unif)

# Alternative models HR (=hazard rate), or unif[orm] :
wren.hr.poly <- ds(wren_lt, key="hr", adjustment="poly", 
                   convert.units=conversion.factor) # simple polynomial adjustment terms
wren.unif.cos <- ds(wren_lt, key="unif", adjustment="cos", # multiple models fitted, succssively adding [cosine] adjmt terms.
                    convert.units=conversion.factor) # ... 3-adj works and retained, with error note for 4-adj (=> don't panic)
# Model comparison, goodness of fit:
AIC(wren.hn, wren.hr.poly, wren.unif.cos)
gof_ds(wren.hr.poly)
# comparison table:
knitr::kable(summarize_ds_models(wren.hn, wren.hr.poly, wren.unif.cos),digits=3,
             caption="Model comparison table for wren line transect data, Montrave.")

# HazardRate "best".  But now look at detn fns:
par(mfrow=c(1,2)) # ="Set graphic PARameters for ... multiframe rows"?
plot(wren.hr.poly, breaks=cutpoints, main="Hazard rate") # p(det) too perfect then too fast to decline
plot(wren.unif.cos, breaks=cutpoints, main="Uniform cosine")

# MILLER et al. code
data("minke")
head(minke); tail(minke)
# Region.Label  Area Sample.Label Effort distance
# 1        South 84734            1  86.75     0.10
# 2        South 84734            1  86.75     0.22
# 3        South 84734            1  86.75     0.16
# ...
# 97        North 630582           24 130.33     0.10
# 98        North 630582           24 130.33     0.10
# 99        North 630582           25 107.72       NA

summary(minke)
# Region.Label            Area         Sample.Label       Effort          distance     
# Length:99          Min.   : 84734   Min.   : 1.00   Min.   :  0.19   Min.   :0.0000  
# Class :character   1st Qu.: 84734   1st Qu.: 5.50   1st Qu.: 75.63   1st Qu.:0.2125  
# Mode  :character   Median :630582   Median :14.00   Median :100.93   Median :0.4900  
#                    Mean   :382469   Mean   :12.73   Mean   : 87.85   Mean   :0.5464  
#                    3rd Qu.:630582   3rd Qu.:18.00   3rd Qu.:126.32   3rd Qu.:0.7950  
#                    Max.   :630582   Max.   :25.00   Max.   :144.90   Max.   :1.9400  
#                                                                      NA's   :9   

minke_hn <- ds(minke, truncation = 1.5) # fits just HN key fn; (tried cosine adjt too)
# 10 iterations = 7.2s; vs 3.87s with no adj: ds(minke, truncation = 1.5, adjustment = NULL) 

plot(minke_hn)
summary(minke_hn)
summinke <- summary(minke_hn)

densum <- summinke$dht$individuals$D # density output, viz.:
#   Label   Estimate          se        cv         lcl        ucl       df
# 1 North 0.02097339 0.007876453 0.3755450 0.009523884 0.04618738 12.27398
# 2 South 0.04681073 0.011281913 0.2410113 0.028272077 0.07750560 15.80275
# 3 Total 0.02403400 0.007179465 0.2987212 0.012838347 0.04499280 14.00459
#= (Stratum)  ^den^      s.e.        c.v.    Lower CI     Upper CI   d.f.


# Fin