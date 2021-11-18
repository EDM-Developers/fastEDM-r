
library(fastEDM)

logistic_map <- function(obs) {
  
  r_x <- 3.625
  r_y <- 3.77
  beta_xy <- 0.05
  beta_yx <- 0.4
  tau <- 1
  
  x <- rep(NA, obs)
  y <- rep(NA, obs)
  
  x[1] <- 0.2
  y[1] <- 0.4
  
  for (i in 2:obs) {
    x[i] <- x[i-1] * (r_x * (1 - x[i-1]) - beta_xy * y[i-1])
    y[i] <- y[i-1] * (r_y * (1 - y[i-1]) - beta_yx * x[i-tau])
  }
  
  return(list(x=x, y=y))
}

obs <- 500
map <- logistic_map(obs)
x <- map$x[300:obs]
y <- map$y[300:obs]
t <- seq_along(x)

# explore x, e(2/10)
cat("Command: explore x, e(2/10)\n\n")
print(edm(t, x, E=2:10)) 

# edm xmap x y, k(5)
cat("Command: edm xmap x y, k(5)\n\n")
print(edm(t, x, y, k=5))

#edm xmap x y, e(6) lib(8)
cat("Command: edm xmap x y, e(6) lib(8)\n\n")
print(edm(t, x, y, E=6, library=8))

#edm explore x, k(5) crossfold(10)
cat("Command: edm explore x, k(5) crossfold(10)\n\n")
print(edm(t, x, k=5, crossfold=10))

# edm explore x, theta(0.2(0.1)2.0) algorithm(smap)
cat("Command: edm explore x, theta(0.2(0.1)2.0) algorithm(smap)\n\n")
print(edm(t, x, theta=seq(0.2, 2.0, 0.1), algorithm="smap")) 

# edm xmap x y, theta(0.2) algorithm(smap) savesmap(beta)
cat("Command: edm xmap x y, theta(0.2) algorithm(smap) savesmap(beta)\n\n")
res <- edm(t, x, y, theta=0.2, algorithm="smap", saveSMAPCoeffs=TRUE)
print(res)

beta <- res$coeffs
stopifnot(sum(is.na(beta)) == 0)

# edm xmap y x, predict(x2) direction(oneway)
cat("Command: edm xmap y x, predict(x2) direction(oneway)\n\n")
res <- edm(t, y, x, savePredictions=TRUE)
print(res)

x2 <- res$predictions
stopifnot(sum(is.na(x2)) == 0) 

# edm explore x, copredict(teste) copredictvar(y)
#cat("Command: edm explore x, copredict(teste) copredictvar(y)\n\n")
#res <- edm(t, x,  copredict(teste) copredictvar(y)
#assert teste!=. if (_n>1 & _n < _N)

# edm explore z.x, p(10)
z.x <- (x - mean(x)) / sd(x)
print(edm(t, z.x, p=10))

# edm xmap y x, p(10) direction(oneway)
print(edm(t, y, x, p=10))

#edm xmap y x, p(10) copredict(testx) copredictvar(x2) direction(oneway)
#assert testx != . if _n >= 3 & _n <= _N - 10

# edm xmap y x, p(10) copredict(testx2) copredictvar(z.x2) direction(oneway)
# assert testx2 != . if _n >= 3 & _n <= _N - 10
# 
# edm xmap y x, extra(u1) p(10) copredict(testx3) copredictvar(z.x2) direction(oneway)
# assert testx3 != . if _n >= 3 & _n <= _N - 10

# Check explore / xmap consistency

# edm xmap l.x x, direction(oneway)
cat("Command: edm xmap l.x x, direction(oneway)\n\n")
l.x <- c(NA, head(x, -1))
resXmap <- edm(t, l.x, x)

# edm explore x, full
cat("Command: edm explore x, full\n\n")
resExplore <- edm(t, x, full=TRUE)

stopifnot(resXmap$summary["rho"] == resExplore$summary["rho"])

# edm explore x, e(2 3) theta(0 1)
cat("Command: edm explore x, e(2 3) theta(0 1)\n\n")
print(edm(t, x, E=c(2, 3), theta=c(0, 1)))

# Test missing data