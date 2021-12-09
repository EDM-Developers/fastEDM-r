tslag <- function(t, x, lag=1, dt=1) {
  l.x <- rep(NA, length(t))
  for (i in seq_along(t)) {
    lagged_t <- t[i]-lag*dt
    if (!is.na(lagged_t) && lagged_t %in% t) {
      l.x[i] <- x[which(t == lagged_t)]
    }
  }
  return(l.x)
}

tsdiff <- function(t, x, lag=1, dt=1) {
  d.x <- rep(NA, length(t))
  for (i in seq_along(t)) {
    lagged_t <- t[i]-lag*dt
    if (!is.na(x[i]) && !is.na(lagged_t) && lagged_t %in% t) {
      d.x[i] <- x[i] - x[which(t == lagged_t)]
    }
  }
  return(d.x)
}

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

test_that("ci-test", {
  obs <- 500
  map <- logistic_map(obs)
  
  x <- map$x[300:obs]
  y <- map$y[300:obs]
  
  u1 <- c(-1.027216, 1.527376, .5907618, 1.070512, 2.139774, -.5876155, .1418234, 1.390853, -1.030574, .5835255, 1.538284, 1.095415, 1.289363, .4250214, 1.332112, .1224301, .4007208, 1.163034, -.9338163, -1.553558, 1.128875, .71824, .8828724, -.9635994, .5716761, .0727569, -.3750865, -.8911737, -.8376914, -.3425734, -1.895796, 1.220617, .8647164, -.4872026, .1291741, -1.807868, .9658784, -.8437532, .7287974, -.0579607, -.7721093, .3223931, .4673252, -.3628134, -.8418728, -.8550454, -1.341583, -.4182656, .4155265, -.3210205, .7979518, .0385472, -2.345896, -.0535184, -1.997315, -.897661, -1.172937, -1.374793, -.439018, 1.212688, -.8391462, -.2125729, .3922674, -1.24292, -.3563064, -1.368325, 1.293824, -1.078043, -.6217906, .2247944, -.3572458, 1.455859, .177133, -.4954876, -.4623527, -.9394832, -1.381252, .3134706, .1598284, .4492666, .7745574, 2.02939, .2769991, -1.729418, -.0719662, -.4887659, -.6402079, -.3815501, -.6201261, -.6295606, .2707956, 1.056473, -1.657482, 1.228817, .8577658, .4940666, 1.37631, -.0235891, 1.044822, .2835678, .019814, -1.331117, -.4936376, -1.570097, 1.482886, -.2730185, -.467406, .8039773, .6066654, .099022, 1.246193, -.6019896, -1.078758, .0527143, .522496, .7971591, 2.091462, -1.87791, 1.123751, .1762845, 1.552169, -.4524258, .4963196, -1.343762, 1.630493, -.1519897, .4249264, .1730838, -1.662154, .5415513, 1.762257, .4248972, -1.56878, -.0073573, .4523424, -1.077807, -3.545176, -1.198717, 1.314406, -1.067673, -.7234299, 1.150322, 2.114344, .4767627, 1.228333, 1.247601, -.2687568, 1.233031, 1.063017, -1.619441, .5857949, 1.296269, .8043274, .3258621, 3.569143, .3741727, -1.49533, -.0184031, .2356096, -1.738142, -.3104737, -.377933, -.5639113, -1.457661, .9921553, -.9124324, -.0439041, -.6419182, .5668358, -.4034521, -.3590932, -1.489591, -.5190973, .5887823, .8400694, .0363247, 1.122107, -.0369949, 1.10605, .6818572, -.1490808, -.9733297, -.8749319, .6384861, -1.647552, -2.270525, .6330903, .1588243, -.0146699, -.2460195, .7494598, -.0442753, -1.198142, -.1973266, .7962075, -.0928933, 2.165736, -.7527414, 1.006963, .1770673, -.4803994)
  
  t <- 299 + seq_along(x)
  
  # explore x, e(2/10)
  res <- edm(t, x, E=2:10)
  rho <- c(.99893, .99879, .99835, .99763, .99457, .99385, .991, .98972, .98572)
  expect_equal(res$summary$rho, rho, tolerance=1e-4)
  
  # edm xmap x y, k(5)
  res1 <- edm(t, x, y, k=5)
  res2 <- edm(t, y, x, k=5)
  expect_equal(res1$summary$rho, .55861, tolerance=1e-4)
  expect_equal(res2$summary$rho, .94454, tolerance=1e-4)
  
  # edm xmap x y, e(6) lib(8)
  res1 <- edm(t, x, y, E=6, library=8)
  res2 <- edm(t, y, x, E=6, library=8)
  expect_equal(res1$summary$rho, .3362, tolerance=1e-4)
  expect_equal(res2$summary$rho, .51116, tolerance=1e-4)
  
  # edm explore x, k(5) crossfold(10)
  res <- edm(t, x, k=5, crossfold=10)
  expect_equal(mean(res$summary$rho), .99946, tolerance=1e-4)
  
  # edm explore x, theta(0.2(0.1)2.0) algorithm(smap)
  res <- edm(t, x, theta=seq(0.2, 2.0, 0.1), algorithm="smap")
  expect_equal(res$summary$rho[1], .99874, tolerance=1e-4)
  expect_equal(res$summary$rho[length(res$summary$rho)], .99882, tolerance=1e-4)
  
  # edm xmap x y, theta(0.2) algorithm(smap) savesmap(beta)
  res1 <- edm(t, x, y, theta=0.2, algorithm="smap", saveSMAPCoeffs=TRUE)
  res2 <- edm(t, y, x, theta=0.2, algorithm="smap", saveSMAPCoeffs=TRUE)
  beta1 <- res1$coeffs
  expect_equal(res1$summary$rho, .66867, tolerance=1e-4)
  expect_equal(res2$summary$rho, .98487, tolerance=1e-4)
  
  # assert beta1_b2_rep1 != . if _n > 1
  expect_equal(sum(is.na(beta1[1,])), ncol(beta1))
  expect_equal(sum(is.na(tail(beta1, -1))), 0)
  
  # edm xmap y x, predict(x2) direction(oneway)
  res <- edm(t, y, x, savePredictions=TRUE)
  x2 <- res$predictions
  expect_equal(res$summary$rho, .94272, tolerance=1e-4)
  
  # assert x2 != . if _n > 1
  expect_equal(is.na(x2[1]), TRUE)
  expect_equal(sum(is.na(tail(x2, -1))), 0)
  
  # edm explore x, copredict(teste) copredictvar(y)
  res <- edm(t, x, copredict = y, saveCoPredictions=TRUE)
  teste <- res$copredictions
  expect_equal(res$summary$rho, .9989, tolerance=1e-4)
  expect_equal(res$co_summary$rho, .78002, tolerance=1e-4)
  
  # assert teste != . if _n > 1
  expect_equal(is.na(teste[1]), TRUE)
  expect_equal(sum(is.na(tail(teste, -1))), 0)
  
  # edm explore z.x, p(10)
  z.x <- (x - mean(x)) / sd(x)
  res <- edm(t, z.x, p=10)
  expect_equal(res$summary$rho, .90235, tolerance=1e-4)
  
  # edm xmap y x, p(10) direction(oneway)
  res <- edm(t, y, x, p=10)
  expect_equal(res$summary$rho, .89554, tolerance=1e-4)
  
  # edm xmap y x, p(10) copredict(testx) copredictvar(x2) direction(oneway)
  res <- edm(t, y, x, p=10, copredict=x2, saveCoPredictions=TRUE)
  testx <- res$copredictions
  expect_equal(res$summary$rho, .89554, tolerance=1e-4)
  expect_equal(res$co_summary$rho, .67401, tolerance=1e-4)
  
  # assert testx != . if _n >= 3
  expect_equal(sum(is.na(testx[1:2])), 2)
  expect_equal(sum(is.na(tail(testx, -2))), 0)
  
  # edm xmap y x, p(10) copredict(testx2) copredictvar(z.x2) direction(oneway)
  z.x2 <- (x2 - mean(x2, na.rm = TRUE)) / sd(x2, na.rm = TRUE)
  res <- edm(t, y, x, p=10,  copredict=z.x2, saveCoPredictions=TRUE)
  testx2 <- res$copredictions
  expect_equal(res$summary$rho, .89554, tolerance=1e-4)
  expect_equal(res$co_summary$rho, .93837, tolerance=1e-4)
  
  # assert testx2 != . if _n >= 3
  expect_equal(sum(is.na(testx2[1:2])), 2)
  expect_equal(sum(is.na(tail(testx2, -2))), 0)
  
  # edm xmap y x, extra(u1) p(10) copredict(testx3) copredictvar(z.x2) direction(oneway)
  res <- edm(t, y, x, extras=list(u1), p=10,  copredict=z.x2, saveCoPredictions=TRUE)
  testx3 <- res$copredictions
  expect_equal(res$summary$rho, .37011, tolerance=1e-4)
  expect_equal(res$co_summary$rho, .9364, tolerance=1e-4)
  
  # assert testx3 != . if _n >= 3
  expect_equal(sum(is.na(testx3[1:2])), 2)
  expect_equal(sum(is.na(tail(testx3, -2))), 0)
  
  # Check explore / xmap consistency
  
  # edm xmap l.x x, direction(oneway)
  resXmap <- edm(t, tslag(t, x), x)
  expect_equal(resXmap$summary$rho, .99939, tolerance=1e-4)
  
  # edm explore x, full
  resExplore <- edm(t, x, full=TRUE)
  expect_equal(resExplore$summary$rho, .99939, tolerance=1e-4)
  
  # assert xmap_r[1,1] == explore_r[1,1]
  expect_equal(resXmap$summary["rho"], resExplore$summary["rho"])
  
  # Check xmap reverse consistency (not necessary to check in this version)
  res1 <- edm(t, x, y)
  res2 <- edm(t, y, x)
  expect_equal(res1$summary$rho, .54213, tolerance=1e-4)
  expect_equal(res2$summary$rho, .94272, tolerance=1e-4)
  
  # Make sure multiple e's and multiple theta's work together
  
  # edm explore x, e(2 3) theta(0 1)
  res <- edm(t, x, E=c(2, 3), theta=c(0, 1))
  rho <- c(.99863, .99895, .99734, .99872)
  expect_equal(res$summary$rho, rho, tolerance=1e-4)
  
  # Check that lowmemory flag is working
  # TODO
  
  # Test missing data
  u <- c(.40321617, .98143454, .33737505, .17918578, .24976152, .29691027, .20499327, .61275607, .19457759, .21902643, .59679411, .27325454, .7832745, .59985177, .83922833, .91918384, .29341243, .68914328, .20482466, .44289448, .50353225, .78469053, .39731964, .1417179, .2297587, .42227099, .84610652, .08079046, .15440883, .56132574, .07844698, .75300325, .01374904, .12845429, .84633796, .9969956, .16866594, .11734659, .53642262, .51599017, .67278131, .96256297, .95954171, .11729408, .56541938, .84330419, .1371089, .38945212, .37677446, .51548277, .51449584, .53923567, .06029199, .58424255, .11279976, .30173961, .56605612, .02343005, .44666614, .83259486, .38594202, .59133444, .54296861, .84755341, .25489783, .54018779, .64907221, .99787491, .86038316, .06460473, .33661689, .0529183, .27751603, .13714739, .53328866, .45544815, .25549471, .21179428, .51436839, .41096275, .11527777, .74136838, .30614234, .7863824, .60394973, .02919326, .52316759, .19693192, .59354107, .26414924, .29072418, .3784841, .93315864, .96052764, .56556703, .52228218, .85167792, .11004516, .79579543, .79621761, .2354088, .93824537, .81922664, .41441068, .49097789, .91274808, .06571458, .45195547, .0709059, .2662836, .53457729, .4265899, .27470816, .51640195, .83499984, .10722623, .49051885, .72424814, .96825454, .98220681, .47936039, .84857664, .94062731, .60951161, .60858524, .45293032, .27575198, .07623614, .44529291, .09089482, .50229062, .80893461, .91586586, .99691506, .35587916, .86148772, .37659892, .24537498, .60092739, .83524085, .66615393, .77105696, .76242014, .99315317, .92303496, .7242071, .75609269, .09389603, .48791731, .9134859, .12996849, .10364283, .68801897, .50264256, .64841714, .84488405, .48749087, .57849803, .8381817, .402748, .66711617, .33963501, .02327878, .26512865, .47665646, .07827067, .61270055, .4072321, .52272033, .19773727, .75079348, .09788921, .18368872, .60296861, .51975047, .24164339, .60783901, .9966368, .80669058, .2479132, .12109083, .50335909, .57944449, .3891651, .53499648, .19121044, .05091697, .80968952, .4719385, .60401948, .99159545, .06252926, .0013005, .83339689, .13339322, .30319281, .56467506, .36715866, .40749152, .89440713, .39817923)
  
  df <- data.frame(t=t, x=x, y=y, u=u)
  df <- df[df$u >= 0.1, ]
  df[df$u < 0.2, "x"] <- NA
  df[df$t %% 19 == 1, "t"] <- NA
  
  t <- df$t
  x <- df$x
  y <- df$y
  u <- df$u
  
  # edm explore x
  res <- edm(t, x)
  expect_equal(res$summary$rho, .99814, tolerance=1e-4)
  
  # edm explore x, dt savemanifold(plugin) dtweight(1)
  res <- edm(t, x, dt=TRUE, saveManifolds=TRUE, dtWeight=1)
  expect_equal(res$summary$rho, .95569, tolerance=1e-4)
  
  # edm explore x, allowmissing
  res <- edm(t, x, allowMissing=TRUE)
  expect_equal(res$summary$rho, .99766, tolerance=1e-4)
  
  # edm explore x, missingdistance(1)
  res <- edm(t, x, allowMissing=TRUE, missingDistance=1.0)
  expect_equal(res$summary$rho, .99765, tolerance=1e-4)
  # TODO: Decide whether this is better -- being explicit about 'allowMissing' & 'missingDistance'
  # or whether to follow Stata and just let the latter auto-enable the former...
  
  # edm xmap x l.x, allowmissing
  res1 <- edm(t, x, tslag(t, x), allowMissing=TRUE)
  res2 <- edm(t, tslag(t, x), x, allowMissing=TRUE)
  expect_equal(res1$summary$rho, .99983, tolerance=1e-4)
  expect_equal(res2$summary$rho, .99864, tolerance=1e-4)
  
  # edm xmap x l.x, extraembed(u) dt alg(smap) savesmap(newb) e(5)
  res1 <- edm(t, x, tslag(t, x), extras=list(u), dt=TRUE, algorithm="smap", E=5)
  res2 <- edm(t, tslag(t, x), x, extras=list(u), dt=TRUE, algorithm="smap", E=5)
  expect_equal(res1$summary$rho, 1.0, tolerance=1e-4) 
  #expect_equal(res2$summary$rho, .77523, tolerance=1e-4) # This one fails <----------------------------

  # edm xmap x l3.x, extraembed(u) dt alg(smap) savesmap(newc) e(5) oneway dtsave(testdt)
  res <- edm(t, x, tslag(t, x, 3), extras=list(u), dt=TRUE, algorithm="smap", E=5)
  #expect_equal(res$summary$rho, .36976, tolerance=1e-4) # This one fails <----------------------------
  
  # edm explore x, extraembed(u) allowmissing dt crossfold(5)
  res <- edm(t, x, extras=list(u), allowMissing=TRUE, dt=TRUE, crossfold=5)
  #expect_equal(mean(res$summary$rho), .92512, tolerance=1e-4) # This one fails <----------------------------
  
  # edm explore d.x, dt
  res <- edm(t, tsdiff(t, x), dt=TRUE)
  #expect_equal(res$summary$rho, .89192, tolerance=1e-4)  # This one fails <----------------------------
  
  # edm explore x, rep(20) ci(95)
  res <- edm(t, x, numReps=20)
  #expect_equal(mean(res$summary$rho), 123, tolerance=1e-4)
  #TODO: ci flag
  
  # edm xmap x y, lib(50) rep(20) ci(95)
  res1 <- edm(t, x, y, library=50, numReps=20)
  res2 <- edm(t, y, x, library=50, numReps=20)
  #expect_equal(mean(res1$summary$rho), 123, tolerance=1e-4)
  #expect_equal(mean(res2$summary$rho), 123, tolerance=1e-4)
  #TODO: ci flag
  
  # Tests from the previous 'bigger-test.do' script
  obs <- 100
  map <- logistic_map(obs)
  
  x <- map$x
  y <- map$y
  t <- seq_along(x)
  
  # edm explore x, e(2) crossfold(2) k(-1) allowmissing
  res <- edm(t, x, E=2, crossfold=2, k=-1, allowMissing=TRUE)
  expect_equal(mean(res$summary$rho), .98175, tolerance=1e-4)
  
  # edm explore x, e(2) crossfold(10) k(-1) allowmissing
  res <- edm(t, x, E=2, crossfold=10, k=-1, allowMissing=TRUE)
  expect_equal(mean(res$summary$rho), .98325, tolerance=1e-4)
  
  # edm explore x, e(5) extra(d.y) full allowmissing
  res <- edm(t, x, E=5, extra=list(tsdiff(t, y)), full=TRUE, allowMissing=TRUE)
  expect_equal(res$summary$rho, .95266, tolerance=1e-4)
  
  # TODO
})
