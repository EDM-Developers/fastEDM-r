library(fastEDM)

df <- read.csv(url(
  "https://github.com/EDM-Developers/fastEDM-r/raw/main/vignettes/chicago.csv"
))

crime_causes_temp <- easy_edm("Crime", "Temperature", data=df, verbosity=0)

temp_causes_crime <- easy_edm("Temperature", "Crime", data=df, verbosity=0)
