library(fastEDM)

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

# Set-up some new edm default arguments for this script
formals(edm)$numThreads <- 4

SAVE_INPUTS <- "ci-test.json"
formals(edm)$saveInputs <- SAVE_INPUTS
if (file.exists(SAVE_INPUTS)) {
  file.remove(SAVE_INPUTS)
}

obs <- 500
map <- logistic_map(obs)

x <- map$x[300:obs]
y <- map$y[300:obs]
t <- 299 + seq_along(x)

# Some normal rv's from Stata using the seed '12345678'
u1 <- c(-1.027216, 1.527376, .5907618, 1.070512, 2.139774, -.5876155, .1418234, 1.390853, -1.030574, .5835255, 1.538284, 1.095415, 1.289363, .4250214, 1.332112, .1224301, .4007208, 1.163034, -.9338163, -1.553558, 1.128875, .71824, .8828724, -.9635994, .5716761, .0727569, -.3750865, -.8911737, -.8376914, -.3425734, -1.895796, 1.220617, .8647164, -.4872026, .1291741, -1.807868, .9658784, -.8437532, .7287974, -.0579607, -.7721093, .3223931, .4673252, -.3628134, -.8418728, -.8550454, -1.341583, -.4182656, .4155265, -.3210205, .7979518, .0385472, -2.345896, -.0535184, -1.997315, -.897661, -1.172937, -1.374793, -.439018, 1.212688, -.8391462, -.2125729, .3922674, -1.24292, -.3563064, -1.368325, 1.293824, -1.078043, -.6217906, .2247944, -.3572458, 1.455859, .177133, -.4954876, -.4623527, -.9394832, -1.381252, .3134706, .1598284, .4492666, .7745574, 2.02939, .2769991, -1.729418, -.0719662, -.4887659, -.6402079, -.3815501, -.6201261, -.6295606, .2707956, 1.056473, -1.657482, 1.228817, .8577658, .4940666, 1.37631, -.0235891, 1.044822, .2835678, .019814, -1.331117, -.4936376, -1.570097, 1.482886, -.2730185, -.467406, .8039773, .6066654, .099022, 1.246193, -.6019896, -1.078758, .0527143, .522496, .7971591, 2.091462, -1.87791, 1.123751, .1762845, 1.552169, -.4524258, .4963196, -1.343762, 1.630493, -.1519897, .4249264, .1730838, -1.662154, .5415513, 1.762257, .4248972, -1.56878, -.0073573, .4523424, -1.077807, -3.545176, -1.198717, 1.314406, -1.067673, -.7234299, 1.150322, 2.114344, .4767627, 1.228333, 1.247601, -.2687568, 1.233031, 1.063017, -1.619441, .5857949, 1.296269, .8043274, .3258621, 3.569143, .3741727, -1.49533, -.0184031, .2356096, -1.738142, -.3104737, -.377933, -.5639113, -1.457661, .9921553, -.9124324, -.0439041, -.6419182, .5668358, -.4034521, -.3590932, -1.489591, -.5190973, .5887823, .8400694, .0363247, 1.122107, -.0369949, 1.10605, .6818572, -.1490808, -.9733297, -.8749319, .6384861, -1.647552, -2.270525, .6330903, .1588243, -.0146699, -.2460195, .7494598, -.0442753, -1.198142, -.1973266, .7962075, -.0928933, 2.165736, -.7527414, 1.006963, .1770673, -.4803994)

# explore x, e(2/10)
cat("\n\nCommand: explore x, e(2/10)\n\n")
res <- edm(t, x, E=2:10)
rho <- c(.99893, .99879, .99835, .99763, .99457, .99385, .991, .98972, .98572)

# edm xmap x y, k(5)
cat("\n\nCommand: edm xmap x y, k(5)\n\n")
res1 <- edm(t, x, y, k=5)
res2 <- edm(t, y, x, k=5)

# edm xmap x y, e(6) lib(8)
cat("\n\nCommand: edm xmap x y, e(6) lib(8)\n\n")
res1 <- edm(t, x, y, E=6, library=8)
res2 <- edm(t, y, x, E=6, library=8)

# edm explore x, k(5) crossfold(10)
cat("\n\nCommand: edm explore x, k(5) crossfold(10)\n\n")
res <- edm(t, x, k=5, crossfold=10)

# edm explore x, theta(0.2(0.1)2.0) algorithm(smap)
cat("\n\nCommand: edm explore x, theta(0.2(0.1)2.0) algorithm(smap)\n\n")
res <- edm(t, x, theta=seq(0.2, 2.0, 0.1), algorithm="smap")

# edm xmap x y, theta(0.2) algorithm(smap) savesmap(beta)
cat("\n\nCommand: edm xmap x y, theta(0.2) algorithm(smap) savesmap(beta)\n\n")
res1 <- edm(t, x, y, theta=0.2, algorithm="smap", saveSMAPCoeffs=TRUE)
res2 <- edm(t, y, x, theta=0.2, algorithm="smap", saveSMAPCoeffs=TRUE)
beta1 <- res1$coeffs

# assert beta1_b2_rep1 != . if _n > 1

# edm xmap y x, predict(x2) direction(oneway)
cat("\n\nCommand: edm xmap y x, predict(x2) direction(oneway)\n\n")
res <- edm(t, y, x, savePredictions=TRUE)
x2 <- res$predictions

# assert x2 != . if _n > 1

# edm explore x, copredict(teste) copredictvar(y)
cat("\n\nCommand: edm explore x, copredict(teste) copredictvar(y)\n\n")
res <- edm(t, x, copredict = y, saveCoPredictions=TRUE)
teste <- res$copredictions

# assert teste != . if _n > 1

# edm explore z.x, p(10)
cat("\n\nCommand: edm explore z.x, p(10)\n\n")
z.x <- (x - mean(x)) / sd(x)
res <- edm(t, z.x, p=10)

# edm xmap y x, p(10) direction(oneway)
cat("\n\nCommand: edm xmap y x, p(10) direction(oneway)\n\n")
res <- edm(t, y, x, p=10)

# edm xmap y x, p(10) copredict(testx) copredictvar(x2) direction(oneway)
cat("\n\nCommand: edm xmap y x, p(10) copredict(testx) copredictvar(x2) direction(oneway)\n\n")
res <- edm(t, y, x, p=10, copredict=x2, saveCoPredictions=TRUE)
testx <- res$copredictions

# assert testx != . if _n >= 3

# edm xmap y x, p(10) copredict(testx2) copredictvar(z.x2) direction(oneway)
cat("\n\nCommand: edm xmap y x, p(10) copredict(testx2) copredictvar(z.x2) direction(oneway)\n\n")
z.x2 <- (x2 - mean(x2, na.rm = TRUE)) / sd(x2, na.rm = TRUE)
res <- edm(t, y, x, p=10,  copredict=z.x2, saveCoPredictions=TRUE)
testx2 <- res$copredictions

# assert testx2 != . if _n >= 3

# edm xmap y x, extra(u1) p(10) copredict(testx3) copredictvar(z.x2) direction(oneway)
cat("\n\nCommand: edm xmap y x, extra(u1) p(10) copredict(testx3) copredictvar(z.x2) direction(oneway)\n\n")
res <- edm(t, y, x, extras=list(u1), p=10,  copredict=z.x2, saveCoPredictions=TRUE)
testx3 <- res$copredictions

# assert testx3 != . if _n >= 3

# Check explore / xmap consistency

# edm xmap l.x x, direction(oneway)
cat("\n\nCommand: edm xmap l.x x, direction(oneway)\n\n")
resXmap <- edm(t, tslag(t, x), x)

# edm explore x, full
cat("\n\nCommand: edm explore x, full\n\n")
resExplore <- edm(t, x, full=TRUE)

# assert xmap_r[1,1] == explore_r[1,1]

# Check xmap reverse consistency (not necessary to check in this version)
res1 <- edm(t, x, y)
res2 <- edm(t, y, x)

# Make sure multiple e's and multiple theta's work together

# edm explore x, e(2 3) theta(0 1)
cat("\n\nCommand: edm explore x, e(2 3) theta(0 1)\n\n")
res <- edm(t, x, E=c(2, 3), theta=c(0, 1))
rho <- c(.99863, .99895, .99734, .99872)

# Check that lowmemory flag is working
cat("\n\nCommand: edm explore x, e(2) lowmemory\n\n")
res <- edm(t, x, lowMemory=TRUE)

# Check that verbosity > 0 is working
res <- edm(t, x, verbosity=1)

# Check that numThreads > 1 is working
res <- edm(t, x, numThreads=4)

# Test missing data
u <- c(.4032161652730771, .9814345444173197, .3373750503685581, .1791857833127429, .2497615187767713, .2969102692610567, .2049932720581052, .6127560679465497, .1945775905414523, .2190264344998417, .5967941090343319, .2732545386445274, .783274501879668, .5998517662013346, .839228334965931, .9191838391294088, .2934124279570279, .6891432790291355, .2048246627394067, .4428944844872612, .5035322526301146, .7846905293507423, .3973196413139125, .1417179025397167, .2297586957409598, .4222709916638957, .8461065233021708, .0807904559242613, .1544088321271263, .5613257424648498, .0784469845635632, .7530032474630725, .0137490402736743, .1284542865245989, .8463379622597907, .9969955966978709, .1686659431226015, .1173465910722692, .5364226240241997, .5159901701260331, .6727813123885479, .9625629665008296, .9595417075442659, .1172940837748137, .5654193787312168, .8433041905223945, .1371089038434833, .3894521245853525, .3767744584921785, .5154827677849309, .5144958419236422, .5392356674361488, .0602919934618357, .5842425522598768, .1127997590422875, .3017396141777563, .566056123747399, .0234300500483061, .446666137907659, .8325948629304355, .3859420229362572, .5913344449911008, .5429686057208226, .8475534057216773, .2548978318715515, .5401877915192624, .6490722076532688, .9978749125996588, .8603831590200534, .0646047336240946, .3366168869966312, .0529183026796373, .2775160253539241, .13714739352711, .5332886602574597, .4554481474468138, .2554947058562992, .2117942810144512, .5143683948513994, .4109627526795842, .1152777732355069, .7413683768433018, .3061423425767477, .7863824045156084, .6039497263555896, .0291932555292082, .5231675914641521, .1969319200248217, .5935410685602686, .2641492447722841, .2907241830474051, .3784841046138104, .9331586394050942, .9605276395327141, .5655670284878619, .5222821795008065, .8516779163207383, .1100451566808257, .7957954273960438, .7962176143105376, .2354088022022441, .9382453662561788, .8192266352790453, .4144106771809849, .490977891655796, .9127480800163904, .0657145839743215, .4519554710752521, .0709059022401829, .2662836024528649, .5345772940845596, .4265899010392383, .2747081587127954, .5164019507468386, .8349998398567205, .107226227911105, .4905188464764981, .7242481387577447, .9682545409150528, .9822068122956905, .4793603880336778, .8485766361764732, .9406273141585527, .6095116145179326, .6085852384563054, .4529303236044436, .2757519773710629, .0762361381946745, .4452929061484255, .0908948151163734, .5022906216805424, .8089346066118283, .9158658613180911, .9969150573511568, .3558791562632609, .8614877244157376, .3765989177730652, .2453749846785717, .6009273882813092, .8352408531141619, .6661539258782491, .7710569642593318, .7624201378973611, .9931531663177428, .923034959330067, .7242071002415561, .7560926864138929, .093896026984115, .4879173119413863, .913485895383184, .1299684863347993, .1036428325657651, .6880189676686744, .5026425620798053, .648417138502822, .8448840509985992, .4874908670167216, .5784980295332124, .8381816958913403, .4027480031037615, .6671161681380354, .3396350096610617, .0232787820504924, .265128654258174, .4766564574754034, .0782706739728657, .6127005528438244, .4072321046516864, .5227203330605766, .1977372747415979, .7507934791611285, .0978892111143869, .1836887170888787, .6029686128691382, .5197504708621133, .2416433912954606, .6078390094509868, .9966368006778291, .8066905771208216, .2479132038263475, .1210908289739182, .5033590914183873, .5794444917633742, .3891651009205163, .5349964783397266, .1912104418458377, .0509169682297838, .8096895160137013, .4719385020310226, .6040194781564009, .9915954480975472, .0625292627821351, .0013004956320272, .8333968865756917, .1333932163423718, .3031928113136991, .5646750635559817, .3671586644732209, .4074915152964259, .8944071321855526, .3981792341304314)

df <- data.frame(t=t, x=x, y=y, u=u)
df <- df[df$u >= 0.1, ]
df[df$u < 0.2, "x"] <- NA
df[df$t %% 19 == 1, "t"] <- NA

t <- df$t
x <- df$x
y <- df$y
u <- df$u

# edm explore x
cat("\n\nCommand: edm explore x\n\n")
res <- edm(t, x)

# edm explore x, dt savemanifold(plugin) dtweight(1)
cat("\n\nCommand: edm explore x, dt savemanifold(plugin) dtweight(1)\n\n")
res <- edm(t, x, dt=TRUE, saveManifolds=TRUE, dtWeight=1)

# edm explore x, allowmissing
cat("\n\nCommand: edm explore x, allowmissing\n\n")
res <- edm(t, x, allowMissing=TRUE)

# edm explore x, missingdistance(1)
cat("\n\nCommand: edm explore x, missingdistance(1)\n\n")
res <- edm(t, x, allowMissing=TRUE, missingDistance=1.0)
# TODO: Decide whether this is better -- being explicit about 'allowMissing' & 'missingDistance'
# or whether to follow Stata and just let the latter auto-enable the former...

# edm xmap x l.x, allowmissing
cat("\n\nCommand: edm xmap x l.x, allowmissing\n\n")
res1 <- edm(t, x, tslag(t, x), allowMissing=TRUE)
res2 <- edm(t, tslag(t, x), x, allowMissing=TRUE)

# edm xmap x l.x, extraembed(u) dt alg(smap) savesmap(newb) e(5)
cat("\n\nCommand: edm xmap x l.x, extraembed(u) dt alg(smap) savesmap(newb) e(5)")
res1 <- edm(t, x, tslag(t, x), extras=list(u), dt=TRUE, algorithm="smap", saveSMAPCoeffs=TRUE, E=5)
res2 <- edm(t, tslag(t, x), x, extras=list(u), dt=TRUE, algorithm="smap", saveSMAPCoeffs=TRUE, E=5)

# edm xmap x l3.x, extraembed(u) dt alg(smap) savesmap(newc) e(5) oneway dtsave(testdt)
cat("\n\nCommand: edm xmap x l3.x, extraembed(u) dt alg(smap) savesmap(newc) e(5) oneway dtsave(testdt)")
res <- edm(t, x, tslag(t, x, 3), extras=list(u), dt=TRUE, algorithm="smap", saveSMAPCoeffs=TRUE, E=5)

# edm explore x, extraembed(u) allowmissing dt crossfold(5)
res <- edm(t, x, extras=list(u), allowMissing=TRUE, dt=TRUE, crossfold=5)
#print(mean(res$summary$rho))

# edm explore d.x, dt
res <- edm(t, tsdiff(t, x), dt=TRUE)

# edm explore x, rep(20) ci(95)
res <- edm(t, x, numReps=20)

# edm xmap x y, lib(50) rep(20) ci(95)
res1 <- edm(t, x, y, library=50, numReps=20)
res2 <- edm(t, y, x, library=50, numReps=20)

# Tests from the previous 'bigger-test.do' script

obs <- 100
map <- logistic_map(obs)

x <- map$x
y <- map$y
t <- seq_along(x)

# Some uniform rv's from Stata using the seed '12345678'
u <- c(.4032161652730771, .9814345444173197, .3373750503685581, .1791857833127429, .2497615187767713, .2969102692610567, .2049932720581052, .6127560679465497, .1945775905414523, .2190264344998417, .5967941090343319, .2732545386445274, .783274501879668, .5998517662013346, .839228334965931, .9191838391294088, .2934124279570279, .6891432790291355, .2048246627394067, .4428944844872612, .5035322526301146, .7846905293507423, .3973196413139125, .1417179025397167, .2297586957409598, .4222709916638957, .8461065233021708, .0807904559242613, .1544088321271263, .5613257424648498, .0784469845635632, .7530032474630725, .0137490402736743, .1284542865245989, .8463379622597907, .9969955966978709, .1686659431226015, .1173465910722692, .5364226240241997, .5159901701260331, .6727813123885479, .9625629665008296, .9595417075442659, .1172940837748137, .5654193787312168, .8433041905223945, .1371089038434833, .3894521245853525, .3767744584921785, .5154827677849309, .5144958419236422, .5392356674361488, .0602919934618357, .5842425522598768, .1127997590422875, .3017396141777563, .566056123747399, .0234300500483061, .446666137907659, .8325948629304355, .3859420229362572, .5913344449911008, .5429686057208226, .8475534057216773, .2548978318715515, .5401877915192624, .6490722076532688, .9978749125996588, .8603831590200534, .0646047336240946, .3366168869966312, .0529183026796373, .2775160253539241, .13714739352711, .5332886602574597, .4554481474468138, .2554947058562992, .2117942810144512, .5143683948513994, .4109627526795842, .1152777732355069, .7413683768433018, .3061423425767477, .7863824045156084, .6039497263555896, .0291932555292082, .5231675914641521, .1969319200248217, .5935410685602686, .2641492447722841, .2907241830474051, .3784841046138104, .9331586394050942, .9605276395327141, .5655670284878619, .5222821795008065, .8516779163207383, .1100451566808257, .7957954273960438, .7962176143105376)

# Some normal rv's from Stata using the seed '1'
u1 <- c(.94173813, .4870331, .55453211, -.57394189, -1.6831859, .20002605, 2.0535631, -1.2874906, .76769561, .57129043, -.9382565, 1.4670297, -2.7969353, .65672988, -.074978352, -.61362195, -1.3412304, .45943514, 1.1464604, 1.3768886, .016770668, .94677925, -.11319048, -.49819016, -1.5304253, -.051611003, -.076513439, -1.3290932, -.45883241, .017877782, .34325397, 1.2092726, .2365011, -.73019648, -.330953, .13359453, 1.0885595, -.63763547, -.42640716, -.014303211, .21588294, .05830165, .059484873, .025059106, 1.0119363, -.35853708, 1.4637038, .70681834, -2.8081942, -.27054599, 1.5580958, .071366407, 2.2807562, .92863506, -.16536251, -.17245923, 2.0830457, -1.6134628, -.16830915, 1.6171873, -.90855205, .0026675737, .82025963, .92624164, 1.6329502, -.232575, -.089815319, -1.0917373, .061252236, 1.1413523, -.0335248, .26932761, -1.9740542, -.99436063, -.53038871, .70026708, -.79605526, -1.1729968, .17358617, -.28859794, .93706262, 1.2917892, -.06885922, 1.0749949, 1.3219627, -.093162067, 1.0999831, .31230453, -.87349302, 1.4867147, -.8970021, -1.1020641, .25990388, -1.9723424, 1.5126398, 1.4318892, -.024286436, -.33137387, -.64844704, -1.7218629)

# edm explore x, e(2) crossfold(2) k(-1) allowmissing
cat("\n\nCommand: edm explore x, e(2) crossfold(2) k(-1) allowmissing\n\n")
res <- edm(t, x, E=2, crossfold=2, k=-1, allowMissing=TRUE)

# edm explore x, e(2) crossfold(10) k(-1) allowmissing
cat("\n\nCommand: edm explore x, e(2) crossfold(10) k(-1) allowmissing\n\n")
res <- edm(t, x, E=2, crossfold=10, k=-1, allowMissing=TRUE)

# edm explore x, e(5) extra(d.y) full allowmissing
cat("\n\nCommand: edm explore x, e(5) extra(d.y) full allowmissing\n\n")
res <- edm(t, x, E=5, extra=list(tsdiff(t, y)), full=TRUE, allowMissing=TRUE)


# Introduce missing data and test all the dt variations
df <- data.frame(t=t, x=x, y=y, u=u, u1=u1)
df <- df[df$u >= 0.1, ]
df[df$u < 0.2, "x"] <- NA
df[df$t %% 7 == 1, "u1"] <- NA
df[df$t %% 19 == 1, "t"] <- NA

t <- df$t
x <- df$x
y <- df$y
u <- df$u
u1 <- df$u1

# Make sure multiple library values are respected

# edm xmap x y, allowmissing dt library(10(5)70)
res1 <- edm(t, x, y, allowMissing=TRUE, dt=TRUE, library=seq(10, 70, 5))
res2 <- edm(t, y, x, allowMissing=TRUE, dt=TRUE, library=seq(10, 70, 5))

# See if the negative values of p are allowed

# edm explore x, p(-1)
res <- edm(t, x, p=-1)

# edm xmap x y, p(-1)
res1 <- edm(t, x, y, p=-1)
res2 <- edm(t, y, x, p=-1)

# Try out copredict and copredictvar combinations with multiple reps etc.

# edm explore x, copredictvar(y)
res <- edm(t, x, copredict=y)

# edm explore x, copredictvar(y) full
res <- edm(t, x, copredict=y, full=TRUE)

# edm xmap x y, copredictvar(u1)
res1 <- edm(t, x, y, copredict=u1)
res2 <- edm(t, y, x, copredict=u1)

# Note the E=5, theta=0 predictions are all the exact same value
# so the correlation being '.' is actually correct.

# edm explore x, e(2/5) theta(0 1) copredictvar(y)
res <- edm(t, x, E=seq(2, 5), theta=c(0, 1), copredict=y)

# edm xmap x y, library(5 10 20 40) copredictvar(u1)
res1 <- edm(t, x, y, library=c(5, 10, 20, 40), copredict=u1)
res2 <- edm(t, y, x, library=c(5, 10, 20, 40), copredict=u1)

# edm explore x, copredictvar(y) rep(20)
res <- edm(t, x, copredict=y, numReps=20)

# edm xmap x y, library(5 10 20 40) copredictvar(u1) rep(100)
res1 <- edm(t, x, y, library=c(5, 10, 20, 40), copredict=u1, numReps=100)
res2 <- edm(t, y, x, library=c(5, 10, 20, 40), copredict=u1, numReps=100)

# # edm explore x, copredictvar(y) rep(100) ci(10)
# res <- edm(t, x, copredict=y, numReps=100)
# 
# # edm xmap x y, library(5 10 20 40) copredictvar(u1) rep(4) detail
# res1 <- edm(t, x, y, library=c(5, 10, 20, 40), copredict=u1, numReps=4)
# res2 <- edm(t, y, x, library=c(5, 10, 20, 40), copredict=u1, numReps=4)

# Panel data
obs <- 100
map <- logistic_map(obs)

x <- map$x
y <- map$y
t <- seq_along(x)
panel <- (seq_along(x) > obs / 3) * 1.0

# edm explore x, e(40)
cat("\n\nCommand: edm explore x, e(40)\n\n")
res <- edm(t, x, panel=panel, E=40)

# edm explore x, e(40) allowmissing
cat("\n\nCommand: edm explore x, e(40) allowmissing\n\n")
res <- edm(t, x, panel=panel, E=40, allowMissing=TRUE)

# edm explore x, e(40) idw(-1)
cat("\n\nCommand: edm explore x, e(40) idw(-1)\n\n")
res <- edm(t, x, panel=panel, E=40, panelWeight=-1)

# edm explore x, e(40) idw(-1) allowmissing
cat("\n\nCommand: edm explore x, e(40) idw(-1) allowmissing\n\n")
res <- edm(t, x, panel=panel, E=40, panelWeight=-1, allowMissing=TRUE)

# edm xmap x y, e(40)
cat("\n\nCommand: edm xmap x y, e(40)\n\n")
res1 <- edm(t, x, y, panel=panel, E=40)
res2 <- edm(t, y, x, panel=panel, E=40)

# edm xmap x y, e(40) allowmissing
cat("\n\nCommand: edm xmap x y, e(40) allowmissing\n\n")
res1 <- edm(t, x, y, panel=panel, E=40, allowMissing=TRUE)
res2 <- edm(t, y, x, panel=panel, E=40, allowMissing=TRUE)

# edm xmap x y, e(40) idw(-1)
cat("\n\nCommand: edm xmap x y, e(40) idw(-1)\n\n")
res1 <- edm(t, x, y, panel=panel, E=40, panelWeight=-1)
res2 <- edm(t, y, x, panel=panel, E=40, panelWeight=-1)

# edm xmap x y, e(40) idw(-1) allowmissing
cat("\n\nCommand: edm xmap x y, e(40) idw(-1) allowmissing\n\n")
res1 <- edm(t, x, y, panel=panel, E=40, panelWeight=-1, allowMissing=TRUE)
res2 <- edm(t, y, x, panel=panel, E=40, panelWeight=-1, allowMissing=TRUE)


# Drop some rows of the dataset & make sure the plugin can handle this
# (i.e. can it replicate a kind of 'tsfill' hehaviour).

# drop if mod(t,7) == 0 
x <- x[t %% 7 != 0]
panel <- panel[t %% 7 != 0]
t <- t[t %% 7 != 0]

# edm explore x, e(5)
cat("\n\nCommand: edm explore x, e(5)\n\n")
res <- edm(t, x, panel=panel, E=5)

# edm explore x, e(5) allowmissing
cat("\n\nCommand: edm explore x, e(5) allowmissing\n\n")
res <- edm(t, x, panel=panel, E=5, allowMissing=TRUE)

# edm explore x, e(5) idw(-1)
cat("\n\nCommand: edm explore x, e(5) idw(-1)\n\n")
res <- edm(t, x, panel=panel, E=5, panelWeight=-1)

# edm explore x, e(5) idw(-1) allowmissing
cat("\n\nCommand: edm explore x, e(5) idw(-1) allowmissing\n\n")
res <- edm(t, x, panel=panel, E=5, panelWeight=-1, allowMissing=TRUE)

# edm explore x, e(5) idw(-1) k(-1)
cat("\n\nCommand: edm explore x, e(5) idw(-1) k(-1)\n\n")
res <- edm(t, x, panel=panel, E=5, panelWeight=-1, k=-1)

# See if the relative dt flags work

# edm explore x, e(5) reldt
cat("\n\nCommand: edm explore x, e(5) reldt\n\n")
res <- edm(t, x, panel=panel, E=5, reldt=TRUE, dtWeight=1, saveManifolds=TRUE)

# edm explore x, e(5) reldt allowmissing
cat("\n\nCommand: edm explore x, e(5) reldt allowmissing\n\n")
res <- edm(t, x, panel=panel, E=5, reldt=TRUE, allowMissing=TRUE)

# edm explore x, e(5) idw(-1) reldt
cat("\n\nCommand: edm explore x, e(5) idw(-1) reldt\n\n")
res <- edm(t, x, panel=panel, E=5, panelWeight=-1, reldt=TRUE)

# edm explore x, e(5) idw(-1) reldt allowmissing
cat("\n\nCommand: edm explore x, e(5) idw(-1) reldt allowmissing\n\n")
res <- edm(t, x, panel=panel, E=5, panelWeight=-1, reldt=TRUE, allowMissing=TRUE)
