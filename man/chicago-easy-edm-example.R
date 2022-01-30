library(fastEDM)
library(readr)

data <- url("https://raw.githubusercontent.com/EDM-Developers/EDM/master/test/chicago.csv")

chicago <- read_csv(data, col_types = cols(crime = col_double()))
chicago <- head(chicago, 500) # Just to speed up the example

crimeCCMCausesTemp <- easy_edm("crime", "temp", data=chicago, verbosity=0)
tempCCMCausesCrime <- easy_edm("temp", "crime", data=chicago, verbosity=0)
