#
# Martin Boer, Biometris, WUR, Wageningen, The Netherlands
#   
# E-mail: martin.boer@wur.nl  
#

#' @export
summary.MMBsplines = function(obj)
{
  cat("lambda:      ", obj$lambda_opt,
      "\neff. dim.    ", obj$ed,
    "\nsigma2:      ", obj$sigma2, 
    "\nb0:          ", obj$a[1],
    "\nb1:          ", obj$a[2], '\n')
  #  "\ncomp. time:  ", obj$time, ' seconds\n')
}

