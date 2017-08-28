#' Maximum variance weights
#'
#' var.weight produces a set of weights that maximizes the total weighted variance of the distribution of different biomarkers within each subject.
#' @param X n by p maxtrix containing observations of p biomarkers of n subjects.
#' @keywords weighting
#' @export
#' @importFrom "stats" "constrOptim"  
#' @examples
#' library(MASS)
#' # a three biomarkers dataset generated from independent normal(0,1)
#' X = mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#' var.weight(X)

var.weight <-function(X){
  p <- ncol(X)

  initial<-rep(1/p, p-1)
  Ci <- c(rep(c(0,-1),p-1),-1)
  Ui <- matrix(0, ncol = p-1, nrow =  2*(p-1))
  for(j in 1:(p-1)){
    Ui[(2*j-1):(2*j), j] <- c(1, -1)
  }
  Ui <- rbind(Ui, rep(-1, p-1))

  cp <- function(w){
    w<-c(w[1:(p-1)], 1-sum(w))
    t1 <- t(w)%*%t(X)^2
    t2 <- t(w)%*%t(X)
    return(sum(t1 - (t2)^2))
  }

  res <- constrOptim(theta = initial, f = cp, grad = NULL, ui = Ui, ci = Ci, control=list(fnscale=-1))$par
  res=c(res,  1- sum(res))
  eval(eval(substitute(expression( res.len <<- length(res) ))))     # set res.len to be used outside the current environment

  res
}
