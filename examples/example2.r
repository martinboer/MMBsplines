#
# Martin Boer, 
# Biometris WUR, Wageningen, The Netherlands.
# 
rm(list = ls())
library('MMBsplines')
packageDescription('MMBsplines')

# set parameters and simulate data:
nobs = 200
xmin = 0
xmax = 10
set.seed(949030)
sim.fun = function(x) { return(3.0 + 0.1*x + sin(2*pi*x))}
x = runif(nobs, min = xmin, max = xmax)
y = sim.fun(x) + 0.5*rnorm(nobs)

# low penalty, high penalty, REML estimate:
df = data.frame(lambda=c(0.001,100.0,1.0),opt = c(FALSE,FALSE,TRUE))

for (i in 1:nrow(df))
{
  obj = MMBsplines(x, y, xmin, xmax, nseg = 100, lambda = df$lambda[i], optimize = df$opt[i])

  # predictions on a dense grid:
  x0 = seq(xmin, xmax, by=0.01)
  yhat = predict(obj, x0)
  ylin = predict(obj, x0, linear=TRUE)
  ysim = sim.fun(x0)
  summary(obj)
  cat('\n')
  # plot:
  title = paste('lambda = ', round(obj$lambda,3))
  plot(x=x, y=y, col='black', pch=16, cex=0.75,ylim=c(0,8), main = title)
  lines(x = x0, y = ysim, col='blue', lwd=2.5)
  lines(x = x0, y = yhat, col='red', lwd=2.5)
  lines(x = x0, y = ylin, col='green', lwd=2.5)
  legend(0.3,8,c("y","ysim","yhat","ylin"), # puts text in the legend
         pch=c(16,NA,NA,NA),
         lty=c(NA,1,1,1), # gives the legend appropriate symbols (lines)
         lwd=c(NA,2.5,2.5,2.5),
         col=c("black","blue","red","green"))   
}  

