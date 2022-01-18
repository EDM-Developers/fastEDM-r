library(fastEDM)
library(readr)

data <- url("https://raw.githubusercontent.com/EDM-Developers/EDM/master/test/chicago.csv")

chicago <- read_csv(data, col_types = cols(crime = col_double()))
chicago <- head(chicago, 500) # Just to speed up the example

easy_edm("crime", "temp", data=chicago)
