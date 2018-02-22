source('~/git/causalTools/causalTools.R')

# set.seed(1162018)
n_obs <- 2500
X <- rnorm(n_obs, sd=1.5)
X <- seq(-4,4,length.out = n_obs)
X <- X[order(X)]
f1 <- 1.75 * X + 1
# f1 <- 3.75 * X + 1
f2 <- 0.1 * X^4  + 0.25*X^2 - 1
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

make_binary <- function(class_labels){
    classes <- as.vector(unique(class_labels))
    classes <- classes[order(classes)]
    n_unique <- length(classes)
    out_mat <- matrix(NA, nrow=length(class_labels), ncol=n_unique)
    for(c in 1:n_unique){
        out_mat[,c] <- ifelse(class_labels == classes[c], 1, 0)
    }
    as.vector(out_mat)
}

make_binary(class_label)

covmat1 <- mc_sqexp_common(as.matrix(X), c(1.25,1.25,1.25))
covmat2 <- mc_normpoly_common(as.matrix(X),
                             c(25,25,25),
                             c(3,3,3), power = 1)
covmat <- covmat1 + covmat2

system.time(outro <- gp_mcla(covmat, make_binary(class_label), 3, max_iters = 50, verbose=F))
system.time(outro2 <- gp_mcla_fast(covmat, make_binary(class_label), 3, max_iters = 50, verbose=F))
all.equal(outro, outro2)

plot(X, p1, ylim=c(0,1), type='l', col='red')
lines(X, p2, ylim=c(0,1), col='blue')
lines(X, p3, ylim=c(0,1), col='green')
points(X, outro$ps[,1], pch=4, col='red')
points(X, outro$ps[,2], pch=4, col='blue')
points(X, outro$ps[,3], pch=4, col='green')



testing <- mcla_optimize(X, make_binary(class_label), mc_sqexp_common, n_classes = 3, c(1))

plot(X, p1, ylim=c(0,1), type='l', col='red')
lines(X, p2, ylim=c(0,1), col='blue')
lines(X, p3, ylim=c(0,1), col='green')
points(X, testing$ps[,1], pch=4, col='red')
points(X, testing$ps[,2], pch=4, col='blue')
points(X, testing$ps[,3], pch=4, col='green')

