set.seed(89123)
n_obs <- 10001
X <- rnorm(n_obs)
X <- as.matrix(X[order(X)])
cov_mat <- gpbalancer::par_sqexp(X, c(1,1))



cov_mat_eigs <- nystrom(cov_mat, n_pts = 500)
est_covmat <- cov_mat_eigs$vectors %*% diag(as.vector(cov_mat_eigs$values)) %*% t(cov_mat_eigs$vectors)


plot(X, cov_mat[n_obs,], type='l')
lines(X, est_covmat[n_obs,], col='red')

apply(abs(cov_mat-est_covmat), 1, mean)

plot(X, apply(abs(cov_mat-est_covmat), 1, mean), type='l')
