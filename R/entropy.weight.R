#' Maximum entropy weights
#'
#' entropy.weight produces a set of weights that maximizes the total weighted entropy of the distribution of different biomarkers within each subject. 
#' @param X n by p maxtrix containing observations of p biomarkers of n subjects.
#' @param h bandwidth for kernel density estimation.
#' @keywords weighting
#' @export
#' @importFrom "stats" "dnorm" "integrate"
#' @examples
#' library(MASS)
#' # a three biomarkers dataset generated from independent normal(0,1)
#' X = mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#' entropy.weight(X, h=1)

entropy.weight <-function(X, h){
#  h <- 1.5
  p <- ncol(X)

  initial<-rep(1/p, p-1)
  Ci <- c(rep(c(0,-1),p-1),-1)
  Ui <- matrix(0, ncol =  p-1, nrow = 2*(p-1))
  for(j in 1:(p-1)){
    Ui[(2*j-1):(2*j), j] <- c(1, -1)
  }
  Ui <- rbind(Ui, rep(-1, p-1))

  X <- t(X)
  cp <- function(w){
    w <- c(w[1:(p-1)],1-sum(w))
    temp <- c(0, ncol(X))

    for (i in 1:ncol(X)) {

      cd <- function(x) {
        c <- 0
        for(j in 1:nrow(X)){
          c <- c + w[j]*dnorm(x, mean = X[j,i], sd = h)
        }
        return(c*log(c))
      }

      temp[i] <- integrate(cd, range(X[,i])[1] - 3*h, range(X[,i])[2] + 3*h)$value
    }
    return(-sum(temp))
  }

  res <- constrOptim(initial, cp, grad = NULL, ui = Ui, ci = Ci, control=list(fnscale=-1))$par
  res=c(res,  1- sum(res))
  eval(eval(substitute(expression( res.len <<- length(res) ))))     # set res.len to be used outside the current environment

  res
}
