library(fastEDM)
library(readr)

chicagoURL <- url("https://github.com/EDM-Developers/fastEDM/raw/master/vignettes/chicago.csv")
chicago <- read.csv(chicagoURL)

crimeCCMCausesTemp <- easy_edm("Crime", "Temperature", data=chicago, verbosity=0)
tempCCMCausesCrime <- easy_edm("Temperature", "Crime", data=chicago, verbosity=0)
