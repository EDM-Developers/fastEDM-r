
# For the progress bar, need the latest Github version of RcppThread:
#   devtools::install_github("tnagler/RcppThread")

library(fastEDM)

opts <- list(nthreads = 1)

t <- c(1, 2, 3, 4)
x <- c(11, 12, 13, 14)

df <- data.frame(t = t, x = x)

es <- 2
libs <- 4

tau <- 1
p <- 1

cat("starting basic test")
res <- run_command(df, es, libs, opts, tau, p, verbosity=1) #, 0, 0, FALSE, FALSE, FALSE)
cat("basic test finished")

#run_command(DataFrame df, NumericVector es, NumericVector libs, List ropts,
#            int tau, int p, int numReps=0, int crossfold=0,
#            bool full=false, bool dtMode=false, bool allowMissing=false)


library(readr)

chicago <- read_csv("chicago.csv", col_types = cols(crime = col_double()))


opts <- list(saveInputs = "testchicago.json", k = -1, nthreads = 8)

df <- data.frame(t = chicago["t"], x = chicago["temp"], y = chicago["crime"])
colnames(df) <- c("t", "x", "y")

libs <- c(10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550,560,570,580,590,600,610,620,630,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800,810,820,830,840,850,860,870,880,890,900,910,920,930,940,950,960,970,980,990,1000,1020,1040,1060,1080,1100,1120,1140,1160,1180,1200,1220,1240,1260,1280,1300,1320,1340,1360,1380,1400,1420,1440,1460,1480,1500,1520,1540,1560,1580,1600,1620,1640,1660,1680,1700,1720,1740,1760,1780,1800,1820,1840,1860,1880,1900,1920,1940,1960,1980,2000,2050,2100,2150,2200,2250,2300,2350,2400,2450,2500,2550,2600,2650,2700,2750,2800,2850,2900,2950,3000,3050,3100,3150,3200,3250,3300,3350,3400,3450,3500,3550,3600,3650,3700,3750,3800,3850,3900,3950,4000,4050,4100,4150,4200,4250,4300,4350,4365)
#libs <- c(10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
#libs <- c(4050,4100,4150,4200,4250,4300,4350,4365)
#libs <- c(10, 100)

library(tictoc)

tic("chicago")
cat("starting chicago test")
res <- run_command(df, es, libs, opts, tau, p, numReps=10, verbosity=0)
cat("finished chicago test")
toc()

# 4 threads 2 tasks 

# 4 threads 430.22 sec , 426.26 sec
# 8 threads broken

# 1 thread 85 secs
# 2 threads 57 secs
# 4 threads 37 secs
# 8 threads 22 secs


library(jsonlite)
parsed <- fromJSON(res)

library(dplyr)

summaryTable <- bind_rows(parsed$summaryTable)

averaged <- aggregate(summaryTable[, c("mae", "rho")], list(summaryTable$library), mean)
colnames(averaged)[[1]] <- "Library"

plot(averaged$Library, averaged$mae)
plot(averaged$Library, averaged$rho)


