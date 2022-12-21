library(fastEDM)

df <- read.csv(url(
  "https://github.com/EDM-Developers/fastEDM-r/raw/main/vignettes/chicago.csv"
))

crimeCCMCausesTemp <- easy_edm("Crime", "Temperature", data=df, verbosity=0)
tempCCMCausesCrime <- easy_edm("Temperature", "Crime", data=df, verbosity=0)
