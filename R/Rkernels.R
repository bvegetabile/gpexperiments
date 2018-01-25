sqexp_poly <- function(X, theta, noise = 1e-4){
    scale0 <- theta[1]
    sig0 <- theta[2]
    scale1 <- 1
    ls1 <- theta[3]


    K1 <- polykernel(X, sig0, pwr = 1,
                     scale = scale0, noise=noise)
    K2 <- sqexp_common(X, lengthscale = ls1,
                       scale = scale1, noise = noise)

    return(K1 + K2)
}
