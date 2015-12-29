#
# Martin Boer, Biometris, WUR, Wageningen, The Netherlands
#   
# E-mail: martin.boer@wur.nl  
#
library(spam)
library(splines)

# converter from splineDesign to spam:
splineDesign.sparse = function(knots, x, ord = 4, derivs, outer.ok = FALSE) {
  as.spam.dgCMatrix(splineDesign(knots, x, ord, derivs, sparse = TRUE, outer.ok = outer.ok))
}

REMLlogprofile.sparse = function(x,obj) {
  lambda = 10^x
  n = obj$N
  p = obj$p
  q = obj$q
  s = n - p
  log.detQ = obj$log_det_Q
  cholC = update(obj$cholC,obj$UtU + lambda * obj$Q_all)
  a = backsolve.spam(cholC, forwardsolve.spam(cholC,obj$Uty))
  yPy = obj$ssy - sum(a*obj$Uty) 
  sigma2 = yPy/(n-p)
  log.detC = 2.0*as.double(determinant(cholC)$modulus)
  logL = -0.5*(log.detC -log.detQ + - q*log(lambda) + s*log(sigma2)  + s)
  logL
}

# non-sparse:
REMLlogprofile.dense = function(x,obj) {
  lambda = 10^x
  n = obj$N
  p = obj$p
  q = obj$q
  s = n - p
  log.detQ = obj$log_det_Q
  C = as.matrix(obj$UtU + lambda * obj$Q_all)
  cholC = chol(C)
  a = backsolve(cholC, forwardsolve(t(cholC),obj$Uty))
  yPy = obj$ssy - sum(a*obj$Uty) 
  sigma2 = yPy/(n-p)
  log.detC = 2.0*as.double(determinant(cholC)$modulus)
  logL = -0.5*(log.detC -log.detQ + - q*log(lambda) + s*log(sigma2)  + s)
  logL
}

#' MMBsplines
#' 
#' @param response name of response variable
#' @param explanatory name of explanatory variable
#' @param data data frame with the data
#' @param xmin minimum value of x
#' @param xmax maximum value of x
#' @param nseg number of segments
#' @param deg degree of B-splines, default deg=2
#' @param sparse sparse model
#' @param lambda penalty parameter
#' @param optimize boolean, default TRUE
#' @param Psplines boolean, default TRUE
#' @export
MMBsplines = function(response, explanatory, data, 
                      xmin, xmax, nseg, deg=2, sparse=TRUE, lambda = 1.0, optimize=TRUE, Psplines=TRUE)
{
  y = data[[response]]
  x = data[[explanatory]]
  t0 = proc.time()[1]
  p = 2
  N = length(y)
  dx = (xmax - xmin) / nseg
  knots = seq(xmin - deg * dx, xmax + deg * dx, by = dx)
  B = splineDesign.sparse(knots, x, derivs=rep(0,N), ord = deg + 1)
  
  m = ncol(B)
  X = matrix(outer(x,c(0:(p-1)),"^"),ncol=p)
  D = diff(diag.spam(m), diff=2)
  if (sparse==FALSE) {
    Z = B %*% t(D) %*% solve(D%*%t(D))
  } else
    Z = B %*% t(D)
  
  U = cbind(X,Z)
  UtU = t(U) %*% U
  Uty = t(U) %*% y
  if (sparse == FALSE)
  {
    Q = diag(m-2)
    log_det_Q = 0.0
  }
  else
  {
    A = diag.spam(m-2)
    if (deg == 3 && !Psplines) {
      A = (1/6)*(abs(row(A)-col(A))==1) + diag(4/6,m-2)
    }
    #P = t(D) %*% A %*% D
    B = D %*% t(D)
    Q = B %*% A %*% B
    log_det_B = as.double(determinant(B)$modulus)
    log_det_A = as.double(determinant(A)$modulus)
    log_det_Q = 2*log_det_B + log_det_A
    #if (nseg > 500) {
    #  Q = Q + 1.0e-12*diag.spam(m-2)
    #}  
  }
  
  #log_det_Q = 2.0*as.double(determinant(D %*%t(D))$modulus)
  Q_all = bdiag.spam(diag.spam(0,2),Q)
  penaltylog10 = seq(-6,6, by=0.1)
  q = ncol(U)-2
  ssy = sum(y^2)
  
  # in the for-loop we use update of Cholesky:
  C = UtU + lambda * Q_all
  cholC = chol(C)
  L = list(C, cholC=cholC,UtU=UtU,Uty=Uty,Q_all = Q_all, 
              ssy=ssy, p=p, q=q, N=N, log_det_Q = log_det_Q,xmin=xmin,xmax=xmax,nseg=nseg,
           deg=deg,knots=knots)
  if (optimize)
  {
    if (sparse==FALSE)
    {
      result = optimize(REMLlogprofile.dense,c(-6,6),tol=1.0e-8,obj=L,maximum=TRUE)
    } else
    {
      result = optimize(REMLlogprofile.sparse,c(-6,6),tol=1.0e-8,obj=L,maximum=TRUE)
    }
  
    lambda_opt = 10^(result$maximum)
    logL_opt = result$objective
  } else {
    if (sparse==FALSE)
    {
      logL_opt = REMLlogprofile.dense(log10(lambda),L)
    } else
    {
      logL_opt = REMLlogprofile.sparse(log10(lambda),L)
    }
    lambda_opt = lambda
  }
  if (sparse==FALSE) {
    C = UtU + lambda_opt * Q_all
    cholC = chol(C)
    a = backsolve(cholC, forwardsolve(t(cholC),Uty))
  }
  else {
    cholC = update(cholC,UtU + lambda_opt * Q_all)
    a = backsolve.spam(cholC, forwardsolve.spam(cholC,Uty))
  }
  yPy = ssy - sum(a*Uty) 
  sigma2 = yPy/(N-p)
  
  t1 = proc.time()[1]
  L$time = as.double(t1-t0)
  L$lambda_opt = lambda_opt
  L$logL_opt = logL_opt
  L$cholC = cholC
  L$a = a
  L$sigma2 = sigma2  
  L$method = ifelse(sparse,"SPARSE","DENSE")
  class(L) = "MMBsplines"
  L
}

logLprofile = function(obj,grid)
{
  sapply(grid,REMLlogprofile,obj)  
}

#' @export
predict = function(obj, x, linear = FALSE)
{
  if (linear) {
    mu = obj$a[1]
    beta = obj$a[2]
    ylin = mu + beta*x
    return(ylin)    
  }
  derivs = rep(0,length(x))  
  if (obj$method == "SPARSE")
  {
    B = splineDesign.sparse(obj$knots, x, derivs=rep(0,length(x)), ord = obj$deg + 1)  
    m = ncol(B)
    X = matrix(outer(x,c(0:(obj$p-1)),"^"),ncol=obj$p)
    D = diff(diag.spam(m), diff=2)
    Z = B %*% t(D)
    U = cbind(X,Z)
    yhat = U %*% obj$a
    yhat
  } else {
    B = splineDesign(obj$knots, x, derivs=rep(0,length(x)), ord = obj$deg + 1)  
    m = ncol(B)
    X = matrix(outer(x,c(0:(obj$p-1)),"^"),ncol=obj$p)
    D = diff(diag(m), diff=2)
    Z = B %*% t(D) %*% solve(D%*%t(D))
    U = cbind(X,Z)
    yhat = U %*% obj$a
    yhat    
  }
}

#' @export
summary.MMBsplines = function(obj)
{
  cat("lambda:      ",    obj$lambda_opt, 
    "\nsigma2:      ", obj$sigma2, 
    "\nb0:          ", obj$a[1],
    "\nb1:          ", obj$a[2], '\n')
  #  "\ncomp. time:  ", obj$time, ' seconds\n')
}

# # simple function to convert from Matrix to spam:
# convertMatrix2spam = function(M) {
#   # save in triplet format:
#   u = summary(M)
#   # right input for spam:
#   L = list(i=as.vector(u$i),j=as.vector(u$j),values=as.vector(u$x))
#   M2 = spam(L)
#   pad(M2) = dim(M)
#   M2
# }

