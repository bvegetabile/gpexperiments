n <- 250
A <- matrix(sqrt(abs(rnorm(n^2))), n, n)
B <- t(A) %*% A
# eigen(B)
# rbenchmark::benchmark(eigen(B),
#                       nystrom(B),
#                       replications = 10)


all.equal(B, nystrom(B, 10))

X <- nystrom(B, 5)
Y <- nystrom(B, 5)

image(B - nystrom(B, 240))

n_obs <- 500
X <- cbind(rnorm(n_obs), rchisq(n_obs, 4))
cov_mat <- gpbalancer::par_sqexp(X, c(1,1,1))
approx_mat <- nystrom(cov_mat, 50)

rbenchmark::benchmark(hm1 <- solve(cov_mat),
                      hm2 <- MASS::ginv(approx_mat),
                      replications=5)

mean(hm1 - hm2)
image(hm1-hm2)
all.equal(cov_mat,approx_mat)

