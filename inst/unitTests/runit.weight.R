test.mdw <- function() {

    library("RUnit")
    library("mdw")
    library("MASS")
    RNGkind("Mersenne-Twister", "Inversion")    
    tolerance=1e-1
    # R.Version()$system is needed b/c 32 bit system gives different results from 64 bit system
    if((file.exists("D:/gDrive/3software/_checkReproducibility") | file.exists("~/_checkReproducibility")) & R.Version()$system %in% c("x86_64, mingw32","x86_64, linux-gnu","x86_64, darwin15.6.0")) tolerance=1e-6 
    print(tolerance)

    ##############################
    # maximum entropy weights
    # generate X, h 
    set.seed(1)
    X <- mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    h <- 1
    w.e <- entropy.weight(X,h)
    v.e <- asym.v.e(X,w.e,h)
    # confirm res
    checkEqualsNumeric(w.e, c(0.3518879, 0.3436133, 0.3044988), tolerance=tolerance)
    checkEqualsNumeric(v.e, matrix(c(0.04372809, -0.01454785, -0.01454785, 0.05235954), byrow = T, ncol = 2), tolerance=tolerance)
    ##############################
    # maximum variance weights
    # generate X
    set.seed(1)
    X <- mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    w.v <- var.weight(X)
    v.v <- asym.v.v(X,w.v)
    # confirm res
    checkEqualsNumeric(w.v, c(0.3593430, 0.3451024, 0.2955546), tolerance=tolerance)
    checkEqualsNumeric(v.v, matrix(c(0.06705141, -0.01533853, -0.01533853, 0.10214535), byrow = T, ncol = 2), tolerance=tolerance)
        



}
