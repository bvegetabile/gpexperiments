n <- 5000
A <- matrix(sqrt(abs(rnorm(n^2))), n, n)
B <- t(A) %*% A

rbenchmark::benchmark(solve(B),
                      eigen(B, symmetric = T),
                      c_eigen(B),
                      nystrom_inv(B, 5),
                      replications = 1)

mean(solve(B) - nystrom_inv(B, 1000, 1e-3))



par(mfrow=c(1,2))
image(solve(B)-nystrom_inv(B, 999))
image(solve(B)-nystrom_inv2(B, 999))

rbenchmark::benchmark(nystrom(B, 500),
                      nystrom_parallel(B, 500),
                      replications = 100)

# n_sims <- 10
# testpts <- seq(25, n, 25)
# outmat <- matrix(NA, nrow=length(testpts), ncol=n_sims)
# for(i in 1:length(testpts)){
#     for(j in 1:n_sims){
#         outmat[i,j]<- mean((B - nystrom(B, testpts[i]))**2)
#     }
# }
# plot(testpts/max(testpts), apply(outmat, 1, mean),
#      ylab = "Mean Squared Error",
#      xlab = "Proportion: m/n")
# plot(log(eigen(B)$values))
#
#
#
# #
# # X <- nystrom(B, 5)
# # Y <- nystrom(B, 5)
# #
# # image(B - nystrom(B, 240))
# #
n_obs <- 1000
prop <- 10
X <- cbind(rnorm(n_obs), rchisq(n_obs, 4))
X <- X[order(X[,2]),]
cov_mat <- gpbalancer::par_sqexp(X, c(1,1,1))

all.equal(nystrom_inv(cov_mat,500), solve(cov_mat))
image(nystrom_inv(cov_mat,1000)-solve(cov_mat))
rbenchmark::benchmark(nystrom_inv(cov_mat,1000),
                      nystrom(cov_mat, 1000),
                      solve(cov_mat), replications=1)

# all.equal(solve(cov_mat), c_inv2(cov_mat))
#
# true_eigs$vectors
# est_eigs$vectors
# approx_mat <- nystrom(cov_mat, 50)
# rbenchmark::benchmark(solve(cov_mat), nystrom(cov_mat, 1000), replications = 1)
#
# mean(cov_mat - est_eigs$vectors %*% diag(as.vector(est_eigs$values)) %*% t(est_eigs$vectors))
# par(mfrow=c(1,1))
# image(cov_mat - est_eigs$vectors %*% diag(as.vector(est_eigs$values)) %*% t(est_eigs$vectors))
#
#
# nystrom(cov_mat, 10)$vectors
# image(nystrom_inv(cov_mat, 500) - solve(cov_mat))
#
# rbenchmark::benchmark(nystrom_inv(cov_mat, 500), solve(cov_mat), replications=1)
# rbenchmark::benchmark(c_inv(cov_mat), solve(cov_mat), replications = 10)
#
#
# rbenchmark::benchmark(eigen(cov_mat), c_eigen(cov_mat), replications = 10)
# rbenchmark::benchmark(solve(cov_mat), replications = 10)
#
# #
# # rbenchmark::benchmark(hm1 <- solve(cov_mat),
# #                       hm2 <- MASS::ginv(approx_mat),
# #                       replications=5)
# #
# # mean(hm1 - hm2)
# # image(hm1-hm2)
# # all.equal(cov_mat,approx_mat)
#
#
# n_sims <- 10
# testpts <- seq(25, n, 25)
# outmat <- matrix(NA, nrow=length(testpts), ncol=n_sims)
# for(i in 1:length(testpts)){
#     for(j in 1:n_sims){
#         outmat[i,j]<- mean((cov_mat - nystrom(cov_mat, testpts[i]))**2)
#     }
# }
# plot(testpts/max(testpts), apply(outmat, 1, mean),
#      ylab = "Mean Squared Error",
#      xlab = "Proportion: m/n")
# plot(eigen(cov_mat)$values)

