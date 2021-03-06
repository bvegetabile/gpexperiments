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


lm_ps<- function(Y, X, wts, true_val = NULL){
    W <- diag(wts)
    invXtWX <- solve(t(X) %*% W %*% X)
    betas <- invXtWX %*% t(X) %*% W %*% Y

    Yhat <- X %*% betas
    resids <- Y - Yhat
    # sighat <- as.double(sum(resids^2) / (length(Y) - 1))
    sighat <- as.double(sum(wts*resids^2) / sum(wts))

    varmat <- sighat * invXtWX %*% t(X) %*% W %*% W %*% X %*% invXtWX

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


lo <- 100

test_rho <- seq(0.05,5,length.out = lo)
n_sims <- 100
n_obs <- 250
true_ate <- 3
mean_results <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
vars_results <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
var2_results <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
bals_results <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
ess1_results <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
ess2_results <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
ess_true <- rep(NA, n_sims)

for(i in 1:n_sims){
    X <- rnorm(n_obs)

    p <- 0.8*pnorm(-0.5 * X^2 + 0.25 * X) + 0.25
    TA <- rbinom(n_obs, 1, p)
    design_mat <- cbind(1, TA)

    wts_true <- make_wts(TA, p)
    wts_t <- wts_true[TA==1]
    wts_c <- wts_true[TA==0]

    ess_true[i] <- sum(wts_t)^2 / sum(wts_t^2) + sum(wts_c)^2 / sum(wts_c^2)
    po_stdev <- 0.25
    YT <- X^2 + 2 + rnorm(n_obs, 0, po_stdev)
    YC <- X + rnorm(n_obs, 0, po_stdev)
    YO <- TA * YT + (1-TA)*YC

    for(r in 1:length(test_rho)){
        cov_mat <- gpbalancer::par_sqexp(as.matrix(X), c(1, test_rho[r]))
        est_obgpps <- gpbalancer::gpbal_fixed(TA, cov_mat, verbose=F)

        wts <- make_wts(TA, est_obgpps$ps)
        bal_obgpps <- gpbalancer::bal_stats(X, TA, wts=wts)
        res <- lm_ps_robust(YO, design_mat, wts, true_ate)

        mean_results[i,r] <- res$ests$coef[2]
        vars_results[i,r] <- res$ests$stderrs[2]^2
        var2_results[i,r] <- lm_ps(YO, design_mat, wts, true_ate)$ests$stderrs[2]^2

        bals_results[i,r] <- bal_obgpps[7]
        ess1_results[i,r] <- sum(wts)^2 / sum(wts^2)

        wts_t <- wts[TA==1]
        wts_c <- wts[TA==0]

        ess2_results[i,r] <- sum(wts_t)^2 / sum(wts_t^2) + sum(wts_c)^2 / sum(wts_c^2)

    }
    message(paste(i,':',sep=""), appendLF = F)
}

lowers <- mean_results - 1.96 * sqrt(vars_results)
uppers <- mean_results + 1.96 * sqrt(vars_results)
covers <- apply((lowers < 3) & (uppers > 3),2, mean)
plot(test_rho, covers, ylim=c(0,1), col=rgb(0,0,0.75,0.5), lwd=5, type='l')

lowers <- mean_results - 1.96 * sqrt(n_obs * var2_results / ess2_results)
uppers <- mean_results + 1.96 * sqrt(n_obs * var2_results / ess2_results)
covers <- apply((lowers < 3) & (uppers > 3),2, mean)
lines(test_rho, covers, col=rgb(0,0.75,0.75,0.5), lwd=5)
abline(h=0.95)

plot(test_rho, apply(mean_results, 2, var) / apply(var2_results, 2, mean))
plot(test_rho,1/apply(n_obs / ess2_results, 2, mean))


pdf('./experiments/simulation_figure_po-sd_three.pdf', height=5, width=8)

plot(0, xlim=range(test_rho), ylim=range(abs(bals_results)),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='Mean Balance - Abs(Std Diff.)',
     main='Covariate Balance')
for(i in 1:n_sims){
    lines(test_rho, abs(bals_results[i,]), lwd=1.5, col=rgb(0,0,0,0.15))
}
lines(test_rho, apply(bals_results, 2, mean),
      lwd=4, col=rgb(0.75,0,0,0.5))
plot(0, xlim=range(test_rho), ylim=range(mean_results),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='ATE Estimate',
     main='Treatment Effect Estimates')
for(i in 1:n_sims){
    lines(test_rho, mean_results[i,], lwd=1.5, col=rgb(0,0,0,0.15))
}
lines(test_rho, apply(mean_results, 2, mean), lwd=4, col=rgb(0.75,0,0,0.5))
abline(h=3, lwd=4, col=rgb(0,0.75,0,0.75))
legend('topright', 'True ATE', lwd=4, col=rgb(0,0.75,0,0.75))
plot(0, xlim=range(test_rho), ylim=c(0,max(vars_results)),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='Estimated Variance',
     main='Robust Variance Estimates')
for(i in 1:n_sims){
    lines(test_rho, vars_results[i,], lwd=1.5, col=rgb(0,0,0,0.15))
}
lines(test_rho, apply(vars_results, 2, mean), lwd=4, col=rgb(0.75,0,0,0.75))
lines(test_rho, apply(mean_results, 2, var), lwd=4, col=rgb(0,0,0.75,0.75))
lines(test_rho, apply(mean_results, 2, var) /  apply(vars_results, 2, mean),
      lwd=4, col=rgb(0.75,0,0.75,0.75))

legend('topright', 'Emp. Var.', lwd=4, col=rgb(0,0,0.75,0.75))
plot(0, xlim=range(test_rho), ylim=c(0,1),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='Coverage Proportion',
     main="Coverage")
lines(test_rho, covers, lwd=3, col=rgb(0,0,0,0.5))
abline(h=0.95, lty=3)
abline(h=1, lty=3)

plot(test_rho, apply(mean_results, 2, var) /  apply(vars_results, 2, mean),
     type='l', lwd=4, col=rgb(0.75,0,0,0.5),
     xlab='Inverse-Length Scale',
     ylab='Ratio',
     main='Ratio of Emp Var to Mean of Robust Var. Est')
abline(h=1, lwd=4, col=rgb(0,0,0,0.5))

dev.off()

plot(apply(ess2_results,2,mean)/n_obs, apply(mean_results, 2, var) /  apply(vars_results, 2, mean))
abline(h=1)
#
#
# dset <- data.frame("YO" = yo, "TA" = t_test, "WTS" = wts_old)
# design.ps <- survey::svydesign(ids=~1, weights=~WTS, data=dset)
# glm1 <- survey::svyglm(YO ~ TA, design=design.ps)
# summary(glm1)$coef[2,1:2]
