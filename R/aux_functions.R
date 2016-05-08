#
# Martin Boer, Biometris, WUR, Wageningen, The Netherlands
#   
# E-mail: martin.boer@wur.nl  
#

# 'converter from splineDesign to spam:
#' @export
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

# not working:
#logLprofile = function(obj,grid)
#{
#  sapply(grid,REMLlogprofile,obj)  
#}
