#
# Martin Boer, 
# Biometris WUR, Wageningen, The Netherlands.
# 
rm(list = ls())
library('MMBsplines')
packageDescription('MMBsplines')

# set parameters and simulate data:
nobs = 1000
xmin = 0
xmax = 10
set.seed(949030)
sim.fun = function(x) { return(3.0 + 0.1*x + sin(2*pi*x))}
x = runif(nobs, min = xmin, max = xmax)
y = sim.fun(x) + 0.5*rnorm(nobs)
data = data.frame(x=x, y=y)
data = data[order(data$x), ]
head(data)

# calculate lambda_max:
obj = MMBsplines(response="y", explanatory="x", data = data, xmin, xmax, nseg = 100)

# predictions on a dense grid:
x0 = seq(xmin, xmax, by=0.01)
yhat = predict(obj, x0)
ylin = predict(obj, x0, linear=TRUE)
ysim = sim.fun(x0)
summary(obj)

# plot:
plot(x=x, y=y, col='black', pch=16, cex=0.5,ylim=c(0,8))
lines(x = x0, y = ysim, col='blue', lwd=2.5)
lines(x = x0, y = yhat, col='red', lwd=2.5)
lines(x = x0, y = ylin, col='green', lwd=2.5)
legend(0.3,8,c("y","ysim","yhat","ylin"), # puts text in the legend
  pch=c(16,NA,NA,NA),
  lty=c(NA,1,1,1), # gives the legend appropriate symbols (lines)
  lwd=c(NA,2.5,2.5,2.5),
  col=c("black","blue","red","green")) 


