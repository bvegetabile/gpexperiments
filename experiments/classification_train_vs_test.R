set.seed(90)
# set.seed(2358)
# set.seed(1349815)
n_obs <- 1000
X <- as.matrix(rnorm(n_obs, sd=1.25))

# fake_cov <- gpbalancer::par_sqexp(X, c(1,0.5))
#
# p <- pnorm(as.vector(mvtnorm::rmvnorm(1, sigma = fake_cov)))

make_p <- function(X){
    p <- (X - 3) * (X - 2) * (X + 4) * (X + 2) * X
    p <- 0.8 * pnorm(3 * (p / (max(p) - min(p)))) + 0.1
    return(p)
}

train_p <- make_p(X)
TA <- rbinom(n_obs, 1, train_p)


X_test <- as.matrix(seq(min(X),max(X),length.out = 250))
test_p <- make_p(X_test)

plot(X, train_p, ylim=c(0,1), pch=19, col=rgb(0,0,0,0.25))
lines(X_test, make_p(X_test))


testing <- gpbal_fixed_predict(y = TA, X_train = X, X_test = X_test,
                               gpbalancer::par_sqexp,
                               theta = c(1.25,0.6), verbose = T)
plot(X, train_p, pch=19, col=rgb(0,0,0,0.5), ylim=c(0,1), xlim=c(min(X),max(X)))
points(X, testing$training_ps, pch=4, col=rgb(0,0,0,0.5))
lines(X_test, testing$test_predictions, col='red')
points(X, TA)

test_rho <- seq(0.1, 3, 0.1)
mse_results <- matrix(NA, ncol=3, nrow=length(test_rho))
for(r in 1:length(test_rho)){
    model_fit <- gpbal_fixed_predict(y = TA, X_train = X, X_test = X_test,
                                    gpbalancer::par_sqexp,
                                    theta = c(1.25,test_rho[r]), verbose = F)
    train_mse <- mean((model_fit$training_ps - train_p)^2)
    test_mse <- mean((model_fit$test_predictions - test_p)^2)
    mse_results[r, 1] <- test_rho[r]
    mse_results[r, 2] <- train_mse
    mse_results[r, 3] <- test_mse
    message('.', appendLF = 'F')
}

par(mfrow=c(1,2))
plot(mse_results[,1], log(mse_results[,2]),
     type='l', col=rgb(0,0,0.75,0.5), lwd=3,
     ylim=c(range(log(mse_results[,2:3]))))
lines(mse_results[,1], log(mse_results[,3]), col=rgb(0.75,0,0.75,0.5), lwd=3)

min_mse <- which(mse_results[,3]==min(mse_results[,3]))

best_fit <- gpbal_fixed_predict(y = TA, X_train = X, X_test = X_test,
                               gpbalancer::par_sqexp,
                               theta = c(1,test_rho[min_mse]), verbose = T)
plot(X, train_p, pch=19, col=rgb(0,0,0,0.5), ylim=c(0,1), xlim=c(min(X),max(X)))
points(X, best_fit$training_ps, pch=4, col=rgb(0,0,0,0.5))
lines(X_test, best_fit$test_predictions, col='red')
points(X, TA)
