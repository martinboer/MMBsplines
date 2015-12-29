#
# Martin Boer, Biometris, WUR, Wageningen, The Netherlands
#   
# E-mail: martin.boer@wur.nl  
#

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

