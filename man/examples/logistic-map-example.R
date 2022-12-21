library(fastEDM)

num_obs <- 500
map <- logistic_map(num_obs)
df <- data.frame(list(x = map$x, y = map$y))

x_causes_y <- easy_edm("x", "y", data=df, verbosity=0)
y_causes_x <- easy_edm("y", "x", data=df, verbosity=0)