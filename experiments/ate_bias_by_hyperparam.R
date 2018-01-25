
set.seed(1209)

make_wts <- function(ta, ps){
    wts <- data.frame(t=(ta/ps) / sum(ta/ps),
                      c=((1-ta)/(1-ps)) / sum((1-ta)/(1-ps)))
    wts <- ifelse(ta==1, wts$t, wts$c)
    return(wts)
}


lm_ps_robust <- function(Y, X, wts, true_val = NULL){
    W <- diag(wts)
    invXtWX <- solve(t(X) %*% W %*% X)
    betas <- invXtWX %*% t(X) %*% W %*% Y

    Yhat <- X %*% betas
    resids <- Y - Yhat
    sighat <- as.double(sum(resids^2) / (length(Y) - 1))

    varmat <- invXtWX %*% t(X) %*% W %*% diag(as.vector(resids)^2) %*% W %*% X %*% invXtWX

    std_errs <- sqrt(diag(varmat))

    low_int <- betas - 1.96 * std_errs
    upp_int <- betas + 1.96 * std_errs

    res <- cbind(betas, std_errs, low_int, upp_int)
    colnames(res) <- c('coef', 'stderrs', 'low95', 'upp95')

    if(!is.null(true_val)){
        cover_test <- res[2,3] < true_val & res[2,4] > true_val
        return(list('ests' = data.frame(res),
                    'covers' = cover_test))
    } else{
        return(list('ests' = res))
    }
}


n_sims <- 100
n_obs <- 500
test_rho <- seq(0.05, 1.5, 0.05)
mean_results_train1 <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
mean_results_train2 <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
mean_results_test1 <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
mean_results_test2 <- matrix(NA, nrow=n_sims, ncol=length(test_rho))

best_mean_results <- matrix(NA, nrow=n_sims, ncol=4)
best_spot_results <- matrix(NA, nrow=n_sims, ncol=4)
best_rho_results <- matrix(NA, nrow=n_sims, ncol=4)

vars_results_train1 <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
vars_results_train2 <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
vars_results_test1 <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
vars_results_test2 <- matrix(NA, nrow=n_sims, ncol=length(test_rho))

bals_results_train1 <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
bals_results_train2 <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
bals_results_test1 <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
bals_results_test2 <- matrix(NA, nrow=n_sims, ncol=length(test_rho))


for(s in 1:n_sims){
    message(paste(s, ':', sep=""), appendLF = 'F')
    X_train <- as.matrix(rnorm(n_obs, sd=1))
    X_test <- as.matrix(rnorm(n_obs, sd=1))
    p_train <- pnorm(-0.6 * X_train)
    p_test <- pnorm(-0.6 * X_test)
    TA_train <- rbinom(n_obs, 1, p_train)
    TA_test <- rbinom(n_obs, 1, p_test)

    design_train <- cbind(1, TA_train)
    design_test <- cbind(1, TA_test)

    YT_train <- X_train^2 + 2 + rnorm(n_obs, sd=0.25)
    YC_train <- X_train + rnorm(n_obs, sd=0.25)
    YT_test <- X_test^2 + 2 + rnorm(n_obs, sd=0.25)
    YC_test <- X_test + rnorm(n_obs, sd=0.25)

    YO_train <- TA_train * YT_train + (1 - TA_train) * YC_train
    YO_test <- TA_test * YT_test + (1 - TA_test) * YC_test

    X_both <- as.matrix(c(X_train, X_test))
    p_both <- c(p_train, p_test)
    TA_both <- c(TA_train, TA_test)

    bal_before_train <- gpbalancer::bal_stats(X_train, TA_train)
    bal_before_test <- gpbalancer::bal_stats(X_test, TA_test)

    bal_results <- matrix(NA, ncol=4, nrow=length(test_rho))
    mean_results <- matrix(NA, ncol=4, nrow=length(test_rho))
    vars_results <- matrix(NA, ncol=4, nrow=length(test_rho))
    ll_results <- matrix(NA, ncol=2, nrow=length(test_rho))

    for(r in 1:length(test_rho)){
        model_fit <- gpbal_fixed_predict(y = TA_train,
                                         X_train = X_train,
                                         X_test = X_both,
                                         gpbalancer::par_sqexp,
                                         theta = c(1,test_rho[r]),
                                         verbose = F)

        wts_train1 <- make_wts(TA_train, model_fit$training_ps)
        wts_train2 <- make_wts(TA_train, model_fit$test_predictions[1:n_obs])
        wts_test1 <- make_wts(TA_test, pnorm(model_fit$test_posterior[(n_obs+1):(2*n_obs)]))
        wts_test2 <- make_wts(TA_test, model_fit$test_predictions[(n_obs+1):(2*n_obs)])

        bal_results[r, 1] <- gpbalancer::bal_stats(X_train, TA_train, wts=wts_train1)[7]
        bal_results[r, 2] <- gpbalancer::bal_stats(X_train, TA_train, wts=wts_train2)[7]
        bal_results[r, 3] <- gpbalancer::bal_stats(X_test, TA_test, wts=wts_test1)[7]
        bal_results[r, 4] <- gpbalancer::bal_stats(X_test, TA_test, wts=wts_test2)[7]

        est_train1 <- lm_ps_robust(YO_train, design_train, wts=wts_train1, 3)
        est_train2 <- lm_ps_robust(YO_train, design_train, wts=wts_train2, 3)
        est_test1 <- lm_ps_robust(YO_test, design_test, wts=wts_test1, 3)
        est_test2 <- lm_ps_robust(YO_test, design_test, wts=wts_test2, 3)

        mean_results[r, ] <- c(est_train1$ests$coef[2],est_train1$ests$coef[2],
                               est_test1$ests$coef[2],est_test1$ests$coef[2])
        vars_results[r, ] <- c(est_train1$ests$stderrs[2],est_train1$ests$stderrs[2],
                               est_test1$ests$stderrs[2],est_test1$ests$stderrs[2])


        ll_results[r, 1] <- model_fit$log_Z_ep

        message('.', appendLF = 'F')
    }

    best_bal_train1 <- which(abs(bal_results[,1]) == min(abs(bal_results[,1])))
    best_bal_train2 <- which(abs(bal_results[,2]) == min(abs(bal_results[,2])))
    best_bal_test1 <- which(abs(bal_results[,3]) == min(abs(bal_results[,3])))
    best_bal_test2 <- which(abs(bal_results[,4]) == min(abs(bal_results[,4])))

    mean_results_train1[s, ] <- mean_results[,1]
    mean_results_train2[s, ] <- mean_results[,2]
    mean_results_test1[s, ] <- mean_results[,3]
    mean_results_test2[s, ] <- mean_results[,4]

    best_spot_results[s, ] <- c(min(best_bal_train1),
                                min(best_bal_train2),
                                min(best_bal_test1),
                                min(best_bal_test2))

    best_mean_results[s, ] <- c(mean_results[min(best_bal_train1), 1],
                                mean_results[min(best_bal_train2), 2],
                                mean_results[min(best_bal_test1), 1],
                                mean_results[min(best_bal_test2), 2])

    best_rho_results[s, ] <- c(test_rho[min(best_bal_train1)],
                               test_rho[min(best_bal_train2)],
                               test_rho[min(best_bal_test1)],
                               test_rho[min(best_bal_test2)])

    vars_results_train1[s, ] <- vars_results[,1]
    vars_results_train2[s, ] <- vars_results[,2]
    vars_results_test1[s, ] <- vars_results[,3]
    vars_results_test2[s, ] <- vars_results[,4]

    bals_results_train1[s, ] <- bal_results[,1]
    bals_results_train2[s, ] <- bal_results[,2]
    bals_results_test1[s, ] <- bal_results[,3]
    bals_results_test2[s, ] <- bal_results[,4]
    message('\n', appendLF = 'F')
}
par(mfrow=c(1,1))
plot(0, xlim=c(0,1.5), ylim=range(abs(bals_results_test1)), pch=19, col=rgb(0,0,0,0))
for(i in 1:n_sims){
    lines(test_rho, abs(bals_results_test1[i,]), lty=1, col=rgb(0,0,0,0.5))
}

plot_res <- function(res, xlabs='', ylabs=''){
    plot(0, xlim=c(0,1.5), ylim=range(res),
         pch=19, col=rgb(0,0,0,0), xlab = xlabs, ylab=ylabs)
    for(i in 1:n_sims){
        lines(test_rho, res[i,], lty=1, col=rgb(0,0,0,0.5))
    }
    abline(h=3, lwd=3, col=rgb(0.75,0,0,0.5))
}
par(mfrow=c(2,2))
plot_res(mean_results_train1)
plot_res(mean_results_train2)
plot_res(mean_results_test1)
plot_res(mean_results_test2)

par(mfrow=c(2,2))
plot_res(vars_results_train1)
plot_res(vars_results_train2)
plot_res(vars_results_test1)
plot_res(vars_results_test2)

n_sims <- 100
compare_sim <- matrix(NA, nrow=n_sims, ncol = 5)
for(i in 1:n_sims){
    message(paste(i, ':', sep=""), appendLF = 'F')
    X_train <- as.matrix(rnorm(n_obs, sd=1))
    p_train <- pnorm(0.6 * X_train)
    TA_train <- rbinom(n_obs, 1, p_train)

    design_train <- cbind(1, TA_train)

    YT_train <- X_train^2 + 2 + rnorm(n_obs, sd=0.25)
    YC_train <- X_train + rnorm(n_obs, sd=0.25)

    YO_train <- TA_train * YT_train + (1 - TA_train) * YC_train

    est_ps <- gpbalancer::gpbal(X_train, TA_train, gpbalancer::par_sqexp,
                                c(1), return_theta = T)
    compare_wts <- make_wts(TA_train, est_ps$ps)
    compare_wts2 <- make_wts(TA_train, p_train)

    compare_res <- lm_ps_robust(YO_train, design_train, wts = compare_wts, 3)
    compare_res2 <- lm_ps_robust(YO_train, design_train, wts = compare_wts2, 3)
    compare_sim[i, 1] <- est_ps$thetas
    compare_sim[i, 2] <- compare_res$ests$coef[2]
    compare_sim[i, 3] <- compare_res$ests$stderrs[2]
    compare_sim[i, 4] <- compare_res2$ests$coef[2]
    compare_sim[i, 5] <- compare_res2$ests$stderrs[2]
}

var(compare_sim[,2]) / mean(compare_sim[,3]^2)
var(compare_sim[,4]) / mean(compare_sim[,5]^2)

apply(compare_sim, 2, mean)



# #-------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------
# X_train <- as.matrix(rnorm(n_obs, sd=1.25))
# X_test <- as.matrix(rnorm(n_obs, sd=1.25))
# p_train <- pnorm(-0.5 * X_train)
# p_test <- pnorm(-0.5 * X_test)
# TA_train <- rbinom(n_obs, 1, p_train)
# TA_test <- rbinom(n_obs, 1, p_test)
#
# design_train <- cbind(1, TA_train)
# design_test <- cbind(1, TA_test)
#
# YT_train <- X_train^2 + 2 + rnorm(n_obs, sd=0.25)
# YC_train <- X_train + rnorm(n_obs, sd=0.25)
# YT_test <- X_test^2 + 2 + rnorm(n_obs, sd=0.25)
# YC_test <- X_test + rnorm(n_obs, sd=0.25)
#
# YO_train <- TA_train * YT_train + (1 - TA_train) * YC_train
# YO_test <- TA_test * YT_test + (1 - TA_test) * YC_test
#
# X_both <- as.matrix(c(X_train, X_test))
# p_both <- c(p_train, p_test)
# TA_both <- c(TA_train, TA_test)
#
# bal_before_train <- gpbalancer::bal_stats(X_train, TA_train)
# bal_before_test <- gpbalancer::bal_stats(X_test, TA_test)
#
#
# bal_results <- matrix(NA, ncol=4, nrow=length(test_rho))
# bal_results2 <- matrix(NA, ncol=4, nrow=length(test_rho))
# mean_results <- matrix(NA, ncol=4, nrow=length(test_rho))
# vars_results <- matrix(NA, ncol=4, nrow=length(test_rho))
# ll_results <- matrix(NA, ncol=2, nrow=length(test_rho))
#
# for(r in 1:length(test_rho)){
#     model_fit <- gpbal_fixed_predict(y = TA_train,
#                                      X_train = X_train,
#                                      X_test = X_both,
#                                      gpbalancer::par_sqexp,
#                                      theta = c(1,test_rho[r]),
#                                      verbose = F)
#
#     wts_train1 <- make_wts(TA_train, model_fit$training_ps)
#     wts_train2 <- make_wts(TA_train, model_fit$test_predictions[1:n_obs])
#     wts_test1 <- make_wts(TA_test, pnorm(model_fit$test_posterior[(n_obs+1):(2*n_obs)]))
#     wts_test2 <- make_wts(TA_test, model_fit$test_predictions[(n_obs+1):(2*n_obs)])
#
#     bal_results[r, 1] <- gpbalancer::bal_stats(X_train, TA_train, wts=wts_train1)[7]
#     bal_results[r, 2] <- gpbalancer::bal_stats(X_train, TA_train, wts=wts_train2)[7]
#     bal_results[r, 3] <- gpbalancer::bal_stats(X_test, TA_test, wts=wts_test1)[7]
#     bal_results[r, 4] <- gpbalancer::bal_stats(X_test, TA_test, wts=wts_test2)[7]
#
#     est_train1 <- lm_ps_robust(YO_train, design_train, wts=wts_train1, 3)
#     est_train2 <- lm_ps_robust(YO_train, design_train, wts=wts_train2, 3)
#     est_test1 <- lm_ps_robust(YO_test, design_test, wts=wts_test1, 3)
#     est_test2 <- lm_ps_robust(YO_test, design_test, wts=wts_test2, 3)
#
#     mean_results[r, ] <- c(est_train1$ests$coef[2],est_train1$ests$coef[2],
#                            est_test1$ests$coef[2],est_test1$ests$coef[2])
#     vars_results[r, ] <- c(est_train1$ests$stderrs[2],est_train1$ests$stderrs[2],
#                            est_test1$ests$stderrs[2],est_test1$ests$stderrs[2])
#
#
#     ll_results[r, 1] <- model_fit$log_Z_ep
#
#     message('.', appendLF = 'F')
# }
# par(mfcol=c(3,2))
# plot(test_rho, abs(bal_results[,1]), type='l', lwd=3, col=rgb(0.75,0,0,0.5))
# lines(test_rho, abs(bal_results[,2]), type='l', lty=3, lwd=3, col=rgb(0.75,0,0,0.5))
#
# plot(test_rho, mean_results[,1], type='l', lwd=3, col=rgb(0.75,0,0,0.5))
# lines(test_rho, mean_results[,2], type='l', lty=3, lwd=3, col=rgb(0.75,0,0,0.5))
# plot(test_rho, vars_results[,1], type='l', lwd=3, col=rgb(0.75,0,0,0.5))
# lines(test_rho, vars_results[,2], type='l', lty=3, lwd=3, col=rgb(0.75,0,0,0.5))
#
# plot(test_rho, abs(bal_results[,3]), type='l', lwd=3, col=rgb(0,0,0.75,0.5))
# lines(test_rho, abs(bal_results[,4]), type='l', lty=3, lwd=3, col=rgb(0,0,0.75,0.5))
# plot(test_rho, mean_results[,3], type='l', lwd=3, col=rgb(0,0,0.75,0.5))
# lines(test_rho, mean_results[,4], type='l', lty=3, lwd=3, col=rgb(0,0,0.75,0.5))
# plot(test_rho, vars_results[,3], type='l', lwd=3, col=rgb(0,0,0.75,0.5))
# lines(test_rho, vars_results[,4], type='l', lty=3, lwd=3, col=rgb(0,0,0.75,0.5))
