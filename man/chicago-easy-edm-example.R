library(fastEDM)
library(readr)

chicagoURL <- url("https://github.com/EDM-Developers/fastEDM/raw/master/vignettes/chicago.csv")
chicago <- read.csv(chicagoURL)

# Reduce the size of the dataset to speed up the example.
# Note, this will invalidate the results just below. 
chicago <- head(chicago, 500) 

crimeCCMCausesTemp <- easy_edm("Crime", "Temperature", data=chicago, verbosity=0)
tempCCMCausesCrime <- easy_edm("Temperature", "Crime", data=chicago, verbosity=0)
