#
# Martin Boer
#
rm(list = ls())
source('MMBsplines.r')

# Set parameters:
nobs = 1000
xmin = 0
xmax = 10
nseg = 100

# Simulate data:
set.seed(949030)
sim.fun = function(x) { return(3.0 + 0.1*x + sin(2*pi*x))}
x = runif(nobs, min=xmin, max=xmax)
y = sim.fun(x) + 0.5*rnorm(nobs)

# run MMBsplines, using quadratic B-splines (default)
obj = MMBsplines(x, y, xmin, xmax, nseg = nseg)

# make predictions on a small grid: 
x0 = seq(xmin, xmax, by=0.01)
yhat = predict(obj, x0)
ylin = predict(obj, x0, linear=TRUE)
ysim = sim.fun(x0)

# print summary:
summary(obj)

# plot:
plot(x=x, y=y, col='black', pch=16, cex=0.5)
lines(x = x0, y = ysim, col='blue', lwd=2.5)
lines(x = x0, y = yhat, col='red', lwd=2.5)
lines(x = x0, y = ylin, col='green', lwd=2.5)

