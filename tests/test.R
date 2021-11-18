
# For the progress bar, need the latest Github version of RcppThread:
#   devtools::install_github("tnagler/RcppThread")

library(fastEDM)
library(tictoc)

cat("Starting a basic test")

t <- c(1, 2, 3, 4)
x <- c(11, 12, 13, 14)

tic()
res <- edm(t, x, E=2, library=4, tau=1, p=1)
toc()

print(res)

cat("Basic test finished\n\n")

cat("Starting the Chicago crime-temp test\n\n")

library(readr)

chicago <- read_csv("chicago.csv", col_types = cols(crime = col_double()))

libs <- c(seq(10, 200, 5), seq(210, 1000, 10), seq(1020, 2000, 20),
          seq(2050, 4350, 50), 4365)

libs <- seq(10, 4365, 5)

tic()
res1 <- edm(chicago["t"], chicago["temp"], chicago["crime"],
                       E=7, library=libs, numReps=1, makePlots=TRUE, verbosity=0, numThreads=4)

res2 <- edm(chicago["t"], chicago["crime"], chicago["temp"],
                       E=7, library=libs, numReps=1, makePlots=TRUE, verbosity=0, numThreads=4)
toc()

cat("Finished chicago test\n\n")


library(ggplot2)

averaged1 <- aggregate(res1$summary[, c("mae", "rho")], list(res1$summary$library), mean)
averaged2 <- aggregate(res2$summary[, c("mae", "rho")], list(res2$summary$library), mean)
colnames(averaged1)[[1]] <- "library"
colnames(averaged2)[[1]] <- "library"

#df <- data.frame(library = averaged1["library"], tempToCrime = averaged1["rho"],
#                 crimeToTemp = averaged2["rho"])

p = ggplot() + 
  geom_line(data = averaged1, aes(x = library, y = rho), color = "blue") +
  geom_line(data = averaged2, aes(x = library, y = rho), color = "red") +
  xlab("Library") +
  ylab("Rho")

print(p)


p = ggplot() + 
  geom_point(data = res1$summary, aes(x = library, y = rho), alpha = 0.05, color = 2) +
  geom_point(data = res2$summary, aes(x = library, y = rho), alpha = 0.05, color = 3) +
  xlab("Library") +
  ylab("Rho") +
  geom_smooth(method = "loess", alpha = 0.9)

print(p)

points(averaged1$library, averaged1$rho, col=1)
plot(averaged2$library, averaged2$rho, col=2)

plot(averaged$library, averaged$rho)
