set.seed(1242018)


make_wts <- function(ta, ps){
    wts <- data.frame(t=(ta/ps) / sum(ta/ps),
                      c=((1-ta)/(1-ps)) / sum((1-ta)/(1-ps)))
    wts <- ifelse(ta==1, wts$t, wts$c)
    return(wts)
}

make_wts2 <- function(ta, ps){
    wts <- data.frame(t=ta / sum(ta),
                      c=((1-ta) * ps /(1-ps)) / sum((1-ta) * ps/(1-ps)))
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

lo <- 30

test_rho <- seq(0.05,2,length.out = lo)
n_sims <- 100
n_obs <- 500
true_ate <- 3

mean_train <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
vars_train <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
bals_train <- matrix(NA, nrow=n_sims, ncol=length(test_rho))

mean_test <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
vars_test <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
bals_test <- matrix(NA, nrow=n_sims, ncol=length(test_rho))


for(i in 1:n_sims){
    X <- as.matrix(rnorm(n_obs))

    p <- 0.8*pnorm(-0.5 * X^2 + 0.25 * X) + 0.25
    TA <- rbinom(n_obs, 1, p)

    sample_split <- sample(1:n_obs, replace = F)
    train_set <- sample_split[1:mid_pt]
    test_set <- sample_split[1:mid_pt]
    mid_pt <- floor(n_obs/2)
    X_train <- X[train_set,]
    T_train <- TA[train_set]

    X_test <- X[test_set,]
    T_test <- TA[test_set]

    design_mat <- cbind(1, TA)

    po_stdev <- 0.25
    YT <- X + 3 + rnorm(n_obs, 0, po_stdev)
    YC <- X + rnorm(n_obs, 0, po_stdev)
    YO <- TA * YT + (1-TA)*YC

    for(r in 1:length(test_rho)){
        est_predict <- gpexperiments::gpbal_fixed_predict(T_train,
                                                          as.matrix(X_train),
                                                          as.matrix(X_test),
                                                          sqexp_common,
                                                          test_rho[r],
                                                          verbose=F)

        wts_train <- make_wts(T_train, est_predict$training_ps)
        wts_test <- make_wts(T_test, est_predict$test_predictions)

        bal_train <- gpbalancer::bal_stats(X_train, T_train, wts=wts_train)[7]
        bal_test <- gpbalancer::bal_stats(X_test, T_test, wts=wts_test)[7]

        ate_train <- lm_ps_robust(YO[train_set], design_mat[train_set,], wts_train, true_ate)
        ate_test <- lm_ps_robust(YO[test_set], design_mat[test_set,], wts_test, true_ate)

        mean_train[i, r] <- ate_train$ests$coef[2]
        vars_train[i, r] <- ate_train$ests$stderrs[2]^2
        bals_train[i, r] <- bal_train

        mean_test[i, r] <- ate_test$ests$coef[2]
        vars_test[i, r] <- ate_test$ests$stderrs[2]^2
        bals_test[i, r] <- bal_test
    }
    message(paste(i,':',sep=""), appendLF = F)
}

l_train <- mean_train - 1.96 * sqrt(vars_train)
u_train <- mean_train + 1.96 * sqrt(vars_train)
c_train <- apply((lowers < 3) & (uppers > 3),2, mean)

l_test <- mean_test - 1.96 * sqrt(vars_test)
u_test <- mean_test + 1.96 * sqrt(vars_test)
c_test <- apply((lowers < 3) & (uppers > 3),2, mean)


par(mfrow=c(1,2))
plot(0, xlim=range(test_rho), ylim=range(abs(bals_train)),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='Mean Balance - Abs(Std Diff.)',
     main='Covariate Balance')
for(i in 1:n_sims){
    lines(test_rho, abs(bals_train[i,]), lwd=1.5, col=rgb(0,0,0,0.15))
}
lines(test_rho, apply(bals_train, 2, mean),
      lwd=4, col=rgb(0.75,0,0,0.5))
plot(0, xlim=range(test_rho), ylim=range(abs(bals_test)),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='Mean Balance - Abs(Std Diff.)',
     main='Covariate Balance')
for(i in 1:n_sims){
    lines(test_rho, abs(bals_test[i,]), lwd=1.5, col=rgb(0,0,0,0.15))
}
lines(test_rho, apply(bals_test, 2, mean),
      lwd=4, col=rgb(0.75,0,0,0.5))

par(mfrow=c(1,2))
plot(0, xlim=range(test_rho), ylim=range(abs(mean_train)),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='Mean Balance - Abs(Std Diff.)',
     main='Estimated ATE')
for(i in 1:n_sims){
    lines(test_rho, abs(mean_train[i,]), lwd=1.5, col=rgb(0,0,0,0.15))
}
lines(test_rho, apply(mean_train, 2, mean),
      lwd=4, col=rgb(0.75,0,0,0.5))
plot(0, xlim=range(test_rho), ylim=range(abs(mean_test)),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='Mean Balance - Abs(Std Diff.)',
     main='Estimated ATE')
for(i in 1:n_sims){
    lines(test_rho, abs(mean_test[i,]), lwd=1.5, col=rgb(0,0,0,0.15))
}
lines(test_rho, apply(mean_test, 2, mean),
      lwd=4, col=rgb(0.75,0,0,0.5))


apply(mean_train,2,var) / apply(vars_train,2, mean)
apply(mean_test,2,var) / apply(vars_test,2, mean)
