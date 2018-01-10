n <- 250
A <- matrix(sqrt(abs(rnorm(n^2))), n, n)
B <- t(A) %*% A

n_sims <- 10
testpts <- seq(25, n, 25)
outmat <- matrix(NA, nrow=length(testpts), ncol=n_sims)
for(i in 1:length(testpts)){
    for(j in 1:n_sims){
        outmat[i,j]<- mean((B - nystrom(B, testpts[i]))**2)
    }
}
plot(testpts/max(testpts), apply(outmat, 1, mean),
     ylab = "Mean Squared Error",
     xlab = "Proportion: m/n")
plot(log(eigen(B)$values))

#
# X <- nystrom(B, 5)
# Y <- nystrom(B, 5)
#
# image(B - nystrom(B, 240))
#
n_obs <- 1000
X <- cbind(rnorm(n_obs), rchisq(n_obs, 4))
cov_mat <- gpbalancer::par_sqexp(X, c(1,1,1))
approx_mat <- nystrom(cov_mat, 50)
#
# rbenchmark::benchmark(hm1 <- solve(cov_mat),
#                       hm2 <- MASS::ginv(approx_mat),
#                       replications=5)
#
# mean(hm1 - hm2)
# image(hm1-hm2)
# all.equal(cov_mat,approx_mat)


n_sims <- 10
testpts <- seq(25, n, 25)
outmat <- matrix(NA, nrow=length(testpts), ncol=n_sims)
for(i in 1:length(testpts)){
    for(j in 1:n_sims){
        outmat[i,j]<- mean((cov_mat - nystrom(cov_mat, testpts[i]))**2)
    }
}
plot(testpts/max(testpts), apply(outmat, 1, mean),
     ylab = "Mean Squared Error",
     xlab = "Proportion: m/n")
plot(eigen(cov_mat)$values)

