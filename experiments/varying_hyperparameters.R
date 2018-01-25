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

lo <- 40

test_rho <- seq(0.05,2,length.out = lo)
test_sig <- seq(0,10, length.out = lo)
n_sims <- 100
n_obs <- 250
true_ate <- 3
mean_results <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
vars_results <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
bals_results <- matrix(NA, nrow=n_sims, ncol=length(test_rho))

poly1_mean_res <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
poly1_vars_res <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
poly1_bals_res <- matrix(NA, nrow=n_sims, ncol=length(test_rho))

poly2_mean_res <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
poly2_vars_res <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
poly2_bals_res <- matrix(NA, nrow=n_sims, ncol=length(test_rho))

poly3_mean_res <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
poly3_vars_res <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
poly3_bals_res <- matrix(NA, nrow=n_sims, ncol=length(test_rho))


for(i in 1:n_sims){
    X <- rnorm(n_obs)

    p <- 0.8*pnorm(-0.5 * X^2 + 0.25 * X) + 0.25
    TA <- rbinom(n_obs, 1, p)
    design_mat <- cbind(1, TA)

    po_stdev <- 0.25
    YT <- X^2 + 2 + rnorm(n_obs, 0, po_stdev)
    YC <- X + rnorm(n_obs, 0, po_stdev)
    YO <- TA * YT + (1-TA)*YC

    for(r in 1:length(test_rho)){
        cov_mat <- gpbalancer::par_sqexp(as.matrix(X), c(1, test_rho[r]))
        est_obgpps <- gpbalancer::gpbal_fixed(TA, cov_mat, verbose=F)

        cov_poly1 <- gpexperiments::polykernel(as.matrix(X), sig_zero = test_sig[r], pwr=1, noise = 1e-6)
        cov_poly2 <- gpexperiments::polykernel(as.matrix(X), sig_zero = test_sig[r], pwr=2, noise = 1e-6)
        cov_poly3 <- gpexperiments::polykernel(as.matrix(X), sig_zero = test_sig[r], pwr=10, noise = 1e-6)

        est_poly1 <- gpbalancer::gpbal_fixed(TA, cov_poly1, verbose=F)
        est_poly2 <- gpbalancer::gpbal_fixed(TA, cov_poly2, verbose=F)
        est_poly3 <- gpbalancer::gpbal_fixed(TA, cov_poly3, verbose=F)

        wts <- make_wts(TA, est_obgpps$ps)
        wts1 <- make_wts(TA, est_poly1$ps)
        wts2 <- make_wts(TA, est_poly2$ps)
        wts3 <- make_wts(TA, est_poly3$ps)


        bal_obgpps <- gpbalancer::bal_stats(X, TA, wts=wts)[7]
        bal_poly1 <- gpbalancer::bal_stats(X, TA, wts=wts1)[7]
        bal_poly2 <- gpbalancer::bal_stats(X, TA, wts=wts2)[7]
        bal_poly3 <- gpbalancer::bal_stats(X, TA, wts=wts3)[7]

        res <- lm_ps_robust(YO, design_mat, wts, true_ate)
        res1 <- lm_ps_robust(YO, design_mat, wts1, true_ate)
        res2 <- lm_ps_robust(YO, design_mat, wts2, true_ate)
        res3 <- lm_ps_robust(YO, design_mat, wts3, true_ate)

        mean_results[i, r] <- res$ests$coef[2]
        vars_results[i, r] <- res$ests$stderrs[2]^2
        bals_results[i, r] <- bal_obgpps

        poly1_mean_res[i, r] <- res1$ests$coef[2]
        poly1_vars_res[i, r] <- res1$ests$stderrs[2]^2
        poly1_bals_res[i, r] <- bal_poly1

        poly2_mean_res[i, r] <- res2$ests$coef[2]
        poly2_vars_res[i, r] <- res2$ests$stderrs[2]^2
        poly2_bals_res[i, r] <- bal_poly2

        poly3_mean_res[i, r] <- res3$ests$coef[2]
        poly3_vars_res[i, r] <- res3$ests$stderrs[2]^2
        poly3_bals_res[i, r] <- bal_poly3

    }
    message(paste(i,':',sep=""), appendLF = F)
}

lowers <- mean_results - 1.96 * sqrt(vars_results)
uppers <- mean_results + 1.96 * sqrt(vars_results)
covers <- apply((lowers < 3) & (uppers > 3),2, mean)

poly_range <- range(c(poly1_mean_res, poly2_mean_res, poly3_mean_res))

plot(test_sig, apply(poly1_mean_res, 2, mean), ylim=poly_range)
lines(test_sig, apply(poly2_mean_res, 2, mean))
lines(test_sig, apply(poly3_mean_res, 2, mean))


plot(test_sig, apply(poly1_mean_res, 2, var) / apply(poly1_vars_res, 2, mean), ylim=c(0.5,1.5))
lines(test_sig, apply(poly2_mean_res, 2, var) / apply(poly2_vars_res, 2, mean))
lines(test_sig, apply(poly3_mean_res, 2, var) / apply(poly3_vars_res, 2, mean))




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

plot(test_rho, apply(mean_results, 2, var) /  apply(vars_results, 2, mean))
abline(h=1)
#
#
# dset <- data.frame("YO" = yo, "TA" = t_test, "WTS" = wts_old)
# design.ps <- survey::svydesign(ids=~1, weights=~WTS, data=dset)
# glm1 <- survey::svyglm(YO ~ TA, design=design.ps)
# summary(glm1)$coef[2,1:2]
