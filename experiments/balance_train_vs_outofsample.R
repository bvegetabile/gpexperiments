# set.seed(2892)
# set.seed(12389)
# set.seed(333)
set.seed(21391023)
# set.seed(25)

# set.seed(14814)

make_wts <- function(ta, ps){
    wts <- data.frame(t=(ta/ps) / sum(ta/ps),
                      c=((1-ta)/(1-ps)) / sum((1-ta)/(1-ps)))
    wts <- ifelse(ta==1, wts$t, wts$c)
    return(wts)
}

n_obs <- 500
X_train <- as.matrix(rnorm(n_obs, sd=1.25))
X_test <- as.matrix(rnorm(n_obs, sd=1.25))
p_train <- pnorm(-0.5 * X_train)
p_test <- pnorm(-0.5 * X_test)
TA_train <- rbinom(n_obs, 1, p_train)
TA_test <- rbinom(n_obs, 1, p_test)

X_both <- as.matrix(c(X_train, X_test))
p_both <- c(p_train, p_test)
TA_both <- c(TA_train, TA_test)

bal_before_both <- gpbalancer::bal_stats(X_both, TA_both)
bal_after_both <- gpbalancer::bal_stats(X_both, TA_both, wts=make_wts(TA_both, p_both))

bal_before_train <- gpbalancer::bal_stats(X_train, TA_train)
bal_before_test <- gpbalancer::bal_stats(X_test, TA_test)

test_rho <- seq(0.1, 1.5, 0.05)
bal_results <- matrix(NA, ncol=5, nrow=length(test_rho))
bal_results2 <- matrix(NA, ncol=4, nrow=length(test_rho))

ll_results <- matrix(NA, ncol=2, nrow=length(test_rho))

for(r in 1:length(test_rho)){
    model_fit <- gpbal_fixed_predict(y = TA_train,
                                     X_train = X_train,
                                     X_test = X_test,
                                     gpbalancer::par_sqexp,
                                     theta = c(1,test_rho[r]),
                                     verbose = F)
    model_fit2 <- gpbal_fixed_predict(y = TA_both,
                                      X_train = X_both,
                                      X_test = X_both,
                                      gpbalancer::par_sqexp,
                                      theta = c(1,test_rho[r]),
                                      verbose = F)

    wts_train <- make_wts(TA_train, model_fit$training_ps)
    wts_test <- make_wts(TA_test, model_fit$test_predictions)
    wts_test2 <- make_wts(TA_test, pnorm(model_fit$test_posterior))

    wts_test3 <- make_wts(TA_both, c(model_fit$training_ps,
                                     model_fit$test_predictions))

    bal_results[r, 1] <- test_rho[r]
    bal_results[r, 2] <- gpbalancer::bal_stats(X_train, TA_train, wts=wts_train)[7]
    bal_results[r, 3] <- gpbalancer::bal_stats(X_test, TA_test, wts=wts_test)[7]
    bal_results[r, 4] <- gpbalancer::bal_stats(X_test, TA_test, wts=wts_test2)[7]
    bal_results[r, 5] <- gpbalancer::bal_stats(X_both, TA_both, wts=wts_test3)[7]

    wts_train <- make_wts(TA_both, model_fit2$training_ps)
    wts_test <- make_wts(TA_both, model_fit2$test_predictions)
    wts_test2 <- make_wts(TA_both, pnorm(model_fit2$test_posterior))

    bal_results2[r, 1] <- test_rho[r]
    bal_results2[r, 2] <- gpbalancer::bal_stats(X_both, TA_both, wts=wts_train)[7]
    bal_results2[r, 3] <- gpbalancer::bal_stats(X_both, TA_both, wts=wts_test)[7]
    bal_results2[r, 4] <- gpbalancer::bal_stats(X_both, TA_both, wts=wts_test2)[7]

    ll_results[r, 1] <- model_fit$log_Z_ep
    ll_results[r, 2] <- model_fit2$log_Z_ep

    message('.', appendLF = 'F')

}

best_bal <- which(abs(bal_results[,3]) == min(abs(bal_results[,3])))

par(mfrow=c(1,3))
plot(bal_results[,1], abs(bal_results[,2]),
     type='l', col=rgb(0,0,0.75,0.5), lwd=3,
     ylim=c(0, max(abs(bal_results[,2:4]))))
lines(bal_results[,1], abs(bal_results[,3]), col=rgb(0.75,0,0,0.5), lwd=3)
lines(bal_results[,1], abs(bal_results[,4]), col=rgb(0,0.75,0.75,0.5), lwd=3)
lines(bal_results[,1],
      abs(bal_results[,2]) + abs(bal_results[,3]),
      col=rgb(0.75,0.75,0.75,0.5), lwd=3)
abline(h=abs(bal_after_both[7]))


plot(bal_results[,1], abs(bal_results[,2]),
     type='l', col=rgb(0,0,0.75,0.5), lwd=3,
     ylim=c(0, max(abs(bal_results[,2:4]))))
lines(bal_results[,1], abs(bal_results[,3]), col=rgb(0.75,0,0,0.5), lwd=3)
lines(bal_results[,1], abs(bal_results[,4]), col=rgb(0,0.75,0.75,0.5), lwd=3)
lines(bal_results[,1],
      (abs(bal_results[,2]) + abs(bal_results[,4]))/2,
      col=rgb(0.15,0.15,0.15,0.75), lwd=3)
lines(bal_results[,1], abs(bal_results[,5]),
      col=rgb(0.8,0.1,0.1,0.75), lwd=3)
abline(h=abs(bal_after_both[7]))


plot(bal_results2[,1], abs(bal_results2[,2]),
     type='l', col=rgb(0,0,0.75,0.5), lwd=3,
     ylim=c(0, max(abs(bal_results[,2:4]))))
lines(bal_results2[,1], abs(bal_results2[,3]), col=rgb(0.75,0,0,0.5), lwd=3)
lines(bal_results2[,1], abs(bal_results2[,4]), col=rgb(0,0.75,0.75,0.5), lwd=3)
abline(h=abs(bal_after_both[7]))
abline(v=test_rho[best_bal])
#
# plot(test_rho, ll_results[,1],
#      type='l', col=rgb(0,0,0.75,0.5), lwd=3,
#      ylim=range(ll_results))
# lines(test_rho, ll_results[,2],
#       type='l', col=rgb(0,0.75,0.75,0.5), lwd=3)
# abline(v=test_rho[best_bal])


