source('~/git/causalTools/causalTools.R')

# set.seed(1162018)
n_obs <- 100
X <- rnorm(n_obs, sd=1.5)
X <- X[order(X)]
f1 <- 1.75 * X + 1
f2 <- 0.1 * X^3  + 0.25*X^2 - 1
p1 <- exp(f1) / (1 + exp(f1) + exp(f2))
p2 <- exp(f2) / (1 + exp(f1) + exp(f2))
p3 <- 1 / (1 + exp(f1) + exp(f2))

plot(X, p1, ylim=c(0,1), type='l', col='red')
lines(X, p2, ylim=c(0,1), col='blue')
lines(X, p3, ylim=c(0,1), col='green')

p <- cbind(p1, p2, p3)

r_sample <- matrix(NA, nrow=n_obs, ncol = 3)
for(i in 1:n_obs){
    r_sample[i, ] <- t(rmultinom(1, size = 1, p[i,]))
}

class_label <- matrix(NA, nrow=n_obs, ncol = 1)
for(j in 1:n_obs){
    class_label[j, ] <- which(r_sample[j,] == 1)
}

covmat <- mc_sqexp_common(as.matrix(X), c(1,1,1))

outro <- gp_mcla(covmat, rep(class_label,3), 3)
