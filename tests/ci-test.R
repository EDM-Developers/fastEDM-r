
suppressMessages(library(fastEDM))

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

u1 <- c(-1.027216, 1.527376, .5907618, 1.070512, 2.139774, -.5876155, .1418234, 1.390853, -1.030574, .5835255, 1.538284, 1.095415, 1.289363, .4250214, 1.332112, .1224301, .4007208, 1.163034, -.9338163, -1.553558, 1.128875, .71824, .8828724, -.9635994, .5716761, .0727569, -.3750865, -.8911737, -.8376914, -.3425734, -1.895796, 1.220617, .8647164, -.4872026, .1291741, -1.807868, .9658784, -.8437532, .7287974, -.0579607, -.7721093, .3223931, .4673252, -.3628134, -.8418728, -.8550454, -1.341583, -.4182656, .4155265, -.3210205, .7979518, .0385472, -2.345896, -.0535184, -1.997315, -.897661, -1.172937, -1.374793, -.439018, 1.212688, -.8391462, -.2125729, .3922674, -1.24292, -.3563064, -1.368325, 1.293824, -1.078043, -.6217906, .2247944, -.3572458, 1.455859, .177133, -.4954876, -.4623527, -.9394832, -1.381252, .3134706, .1598284, .4492666, .7745574, 2.02939, .2769991, -1.729418, -.0719662, -.4887659, -.6402079, -.3815501, -.6201261, -.6295606, .2707956, 1.056473, -1.657482, 1.228817, .8577658, .4940666, 1.37631, -.0235891, 1.044822, .2835678, .019814, -1.331117, -.4936376, -1.570097, 1.482886, -.2730185, -.467406, .8039773, .6066654, .099022, 1.246193, -.6019896, -1.078758, .0527143, .522496, .7971591, 2.091462, -1.87791, 1.123751, .1762845, 1.552169, -.4524258, .4963196, -1.343762, 1.630493, -.1519897, .4249264, .1730838, -1.662154, .5415513, 1.762257, .4248972, -1.56878, -.0073573, .4523424, -1.077807, -3.545176, -1.198717, 1.314406, -1.067673, -.7234299, 1.150322, 2.114344, .4767627, 1.228333, 1.247601, -.2687568, 1.233031, 1.063017, -1.619441, .5857949, 1.296269, .8043274, .3258621, 3.569143, .3741727, -1.49533, -.0184031, .2356096, -1.738142, -.3104737, -.377933, -.5639113, -1.457661, .9921553, -.9124324, -.0439041, -.6419182, .5668358, -.4034521, -.3590932, -1.489591, -.5190973, .5887823, .8400694, .0363247, 1.122107, -.0369949, 1.10605, .6818572, -.1490808, -.9733297, -.8749319, .6384861, -1.647552, -2.270525, .6330903, .1588243, -.0146699, -.2460195, .7494598, -.0442753, -1.198142, -.1973266, .7962075, -.0928933, 2.165736, -.7527414, 1.006963, .1770673, -.4803994)

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
beta <- res$coeffs
print(res)

# assert beta1_b2_rep1 != . if _n > 1
stopifnot(sum(is.na(beta[1,])) == ncol(beta))
stopifnot(sum(is.na(tail(beta, -1))) == 0)

# edm xmap y x, predict(x2) direction(oneway)
cat("Command: edm xmap y x, predict(x2) direction(oneway)\n\n")
res <- edm(t, y, x, savePredictions=TRUE)
x2 <- res$predictions
print(res)

# assert x2 != . if _n > 1
stopifnot(is.na(x2[1]) && sum(is.na(tail(x2, -1))) == 0)

# edm explore x, copredict(teste) copredictvar(y)
cat("Command: edm explore x, copredict(teste) copredictvar(y)\n\n")
res <- edm(t, x, copredict = y, saveCoPredictions=TRUE)
teste <- res$copredictions
print(res)

# assert teste != . if _n > 1
stopifnot(is.na(teste[1]) && sum(is.na(tail(teste, -1))) == 0)

# edm explore z.x, p(10)
z.x <- (x - mean(x)) / sd(x)
print(edm(t, z.x, p=10))

# edm xmap y x, p(10) direction(oneway)
print(edm(t, y, x, p=10))

# edm xmap y x, p(10) copredict(testx) copredictvar(x2) direction(oneway)
res <- edm(t, y, x, p=10, copredict=x2, saveCoPredictions=TRUE)
testx <- res$copredictions
print(res)

stopifnot(sum(is.na(testx[1:2])) == 2)
stopifnot(sum(is.na(tail(testx, -2))) == 0) 

# edm xmap y x, p(10) copredict(testx2) copredictvar(z.x2) direction(oneway)
z.x2 <- (x2 - mean(x2, na.rm = TRUE)) / sd(x2, na.rm = TRUE)
res <- edm(t, y, x, p=10,  copredict=z.x2, saveCoPredictions=TRUE)
testx2 <- res$copredictions
print(res)

# assert testx2 != . if _n >= 3
stopifnot(sum(is.na(testx2[1:2])) == 2)
stopifnot(sum(is.na(tail(testx2, -2))) == 0)

# edm xmap y x, extra(u1) p(10) copredict(testx3) copredictvar(z.x2) direction(oneway)
res <- edm(t, y, x, extras=list(u1), p=10,  copredict=z.x2, saveCoPredictions=TRUE)
testx3 <- res$copredictions

# assert testx3 != . if _n >= 3
stopifnot(sum(is.na(testx3[1:2])) == 2)
stopifnot(sum(is.na(tail(testx3, -2))) == 0)

# Check explore / xmap consistency

# edm xmap l.x x, direction(oneway)
cat("Command: edm xmap l.x x, direction(oneway)\n\n")
l.x <- c(NA, head(x, -1))
resXmap <- edm(t, l.x, x)

# edm explore x, full
cat("Command: edm explore x, full\n\n")
resExplore <- edm(t, x, full=TRUE)

# assert xmap_r[1,1] == explore_r[1,1]
stopifnot(resXmap$summary["rho"] == resExplore$summary["rho"])

# Make sure multiple e's and multiple theta's work together

# edm explore x, e(2 3) theta(0 1)
cat("Command: edm explore x, e(2 3) theta(0 1)\n\n")
print(edm(t, x, E=c(2, 3), theta=c(0, 1)))

# Test missing data

# TODO

# Tests from the previous 'bigger-test.do' script

obs <- 100
map <- logistic_map(obs)

x <- map$x
y <- map$y
t <- seq_along(x)

# edm explore x, e(2) crossfold(2) k(-1) allowmissing
cat("Command: edm explore x, e(2) crossfold(2) k(-1) allowmissing\n\n")
print(edm(t, x, E=2, crossfold=2, k=-1, allowMissing=TRUE))

# edm explore x, e(2) crossfold(10) k(-1) allowmissing
cat("Command: edm explore x, e(2) crossfold(10) k(-1) allowmissing\n\n")
print(edm(t, x, E=2, crossfold=10, k=-1, allowMissing=TRUE))

# edm explore x, e(5) extra(d.y) full allowmissing
cat("Command: edm explore x, e(5) extra(d.y) full allowmissing\n\n")
d.y <- c(NA, diff(y))
print(edm(t, x, E=5, extra=list(d.y), full=TRUE, allowMissing=TRUE))

# Test e-varying extra is the same as specifying the individual lagged extras

# TODO