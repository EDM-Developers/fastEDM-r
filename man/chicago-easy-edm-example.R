library(fastEDM)

chicago <- read.csv(url("https://github.com/EDM-Developers/fastEDM/raw/master/vignettes/chicago.csv"))

crimeCCMCausesTemp <- easy_edm("Crime", "Temperature", data=chicago, verbosity=0)
tempCCMCausesCrime <- easy_edm("Temperature", "Crime", data=chicago, verbosity=0)
