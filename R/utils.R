summary.abslope <- function(object, ...) {
  cat('Adaptive Bayesian slope model\n\n')
  cat('Selected coefficients:\n')
  cat(object$selected, '\n\n')
  cat('Coefficients:\n')
  line <- format(round(object$beta, 4), nsmall=2)
  cat(line)
}

plot.abslope <- function(object, ...) {
  
  n <- nrow(object$betas)
  
  plot(1, type='n',
       xlab='Iteration',
       ylab='Coefficient',
       xlim=c(1, n), ylim=c(-max(abs(object$betas)), max(abs(object$betas))))
  grid()
  
  for(i in object$selected) {
    lines(1:n, object$betas[,i],
          lwd=2,
          col='red')
  }
}
