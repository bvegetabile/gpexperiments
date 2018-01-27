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
true_ate <- .25^2 + 2.5^2 + 2 - 2.5

mean_gp <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
vars_gp <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
var2_gp <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
bals_gp <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
ess_gp <- matrix(NA, nrow=n_sims, ncol=length(test_rho))

mean_glm <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
vars_glm <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
var2_glm <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
bals_glm <- matrix(NA, nrow=n_sims, ncol=length(test_rho))
ess_glm <- matrix(NA, nrow=n_sims, ncol=length(test_rho))


for(i in 1:n_sims){
    n_obs <- 500
    Xt <- rnorm(n_obs/5, mean=2.5, sd=0.25)
    Xc <- rnorm(n_obs - n_obs/5, mean=0, sd=2)
    X <- c(Xt, Xc)
    TA <- c(rep(1, n_obs/5), rep(0, n_obs - n_obs/5))
    design_mat <- cbind(1, TA)

    po_stdev <- 0.25
    YT <- X^2 + 2 + rnorm(n_obs, 0, po_stdev)
    YC <- X + rnorm(n_obs, 0, po_stdev)
    YO <- TA * YT + (1-TA)*YC

    for(r in 1:length(test_rho)){
        cov_mat <- gpbalancer::par_sqexp(as.matrix(X), c(1, test_rho[r]))

        est_obgpps <- gpbalancer::gpbal_fixed(TA, cov_mat, verbose=F)

        wts_gp <- make_wts2(TA, est_obgpps$ps)
        bal_gp <- gpbalancer::bal_stats(X, TA, wts=wts_gp)
        res_gp <- lm_ps_robust(YO, design_mat, wts_gp, true_ate)

        mean_gp[i,r] <- res_gp$ests$coef[2]
        vars_gp[i,r] <- res_gp$ests$stderrs[2]^2
        var2_gp[i,r] <- lm_ps(YO, design_mat, wts_gp, true_ate)$ests$stderrs[2]^2
        bals_gp[i,r] <- bal_gp[7]

        wts_gp_t <- wts_gp[TA==1]
        wts_gp_c <- wts_gp[TA==0]

        ess_gp[i,r] <- sum(wts_gp_t)^2 / sum(wts_gp_t^2) + sum(wts_gp_c)^2 / sum(wts_gp_c^2)

        est_glm <- glm(TA ~ X + I(X^2), family='binomial')

        wts_glm <- make_wts2(TA, fitted(est_glm))
        bal_glm <- gpbalancer::bal_stats(X, TA, wts=wts_glm)
        res_glm <- lm_ps_robust(YO, design_mat, wts_glm, true_ate)

        mean_glm[i,r] <- res_glm$ests$coef[2]
        vars_glm[i,r] <- res_glm$ests$stderrs[2]^2
        var2_glm[i,r] <- lm_ps(YO, design_mat, wts_glm, true_ate)$ests$stderrs[2]^2
        bals_glm[i,r] <- bal_glm[7]

        wts_glm_t <- wts_glm[TA==1]
        wts_glm_c <- wts_glm[TA==0]

        ess_glm[i,r] <- sum(wts_glm_t)^2 / sum(wts_glm_t^2) + sum(wts_glm_c)^2 / sum(wts_glm_c^2)
    }
    message(paste(i,':',sep=""), appendLF = F)
}

lowers <- mean_gp - 1.96 * sqrt(vars_gp)
uppers <- mean_gp + 1.96 * sqrt(vars_gp)
covers <- apply((lowers < true_ate) & (uppers > true_ate),2, mean)
plot(test_rho, covers, ylim=c(0,1), col=rgb(0,0,0.75,0.5), lwd=5, type='l')

lowers1 <- mean_gp - 1.96 * sqrt(ess_gp * vars_gp / n_obs)
uppers1 <- mean_gp + 1.96 * sqrt(ess_gp * vars_gp / n_obs)
covers1 <- apply((lowers1 < true_ate) & (uppers1 > true_ate),2, mean)
lines(test_rho, covers1, ylim=c(0,1), col=rgb(0.75,0,0.75,0.5), lwd=5, type='l')

emp_var <- matrix(apply(mean_gp, 2, var), nrow=n_sims, ncol = length(test_rho), byrow = T)

lowers2 <- mean_gp - 1.96 * sqrt(emp_var)
uppers2 <- mean_gp + 1.96 * sqrt(emp_var)
covers2 <- apply((lowers2 < true_ate) & (uppers2 > true_ate),2, mean)
lines(test_rho, covers2, ylim=c(0,1), col=rgb(0,0.75,0,0.5), lwd=5, type='l')
abline(h=0.95)

lowers3 <- mean_gp - 1.96 * sqrt(var2_gp)
uppers3 <- mean_gp + 1.96 * sqrt(var2_gp)
covers3 <- apply((lowers3 < true_ate) & (uppers3 > true_ate),2, mean)
lines(test_rho, covers3, ylim=c(0,1), col=rgb(0.75,0.75,0.75,0.5), lwd=5, type='l')

lowers4 <- mean_gp - 1.96 * sqrt(ess_gp * var2_gp / n_obs)
uppers4 <- mean_gp + 1.96 * sqrt(ess_gp * var2_gp / n_obs)
covers4 <- apply((lowers4 < true_ate) & (uppers4 > true_ate),2, mean)
lines(test_rho, covers4, ylim=c(0,1), col=rgb(0,0,0,0.75), lwd=5, type='l')

# pdf('./experiments/simulation_figure_po-sd_three.pdf', height=5, width=8)

plot(0, xlim=range(test_rho), ylim=range(abs(bals_gp)),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='Mean Balance - Abs(Std Diff.)',
     main='Covariate Balance')
for(i in 1:n_sims){
    lines(test_rho, abs(bals_gp[i,]), lwd=1.5, col=rgb(0,0,0,0.15))
}
lines(test_rho, apply(bals_gp, 2, mean),
      lwd=4, col=rgb(0.75,0,0,0.5))
plot(0, xlim=range(test_rho), ylim=range(mean_gp),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='ATE Estimate',
     main='Treatment Effect Estimates')
for(i in 1:n_sims){
    lines(test_rho, mean_gp[i,], lwd=1.5, col=rgb(0,0,0,0.15))
}
lines(test_rho, apply(mean_gp, 2, mean), lwd=4, col=rgb(0.75,0,0,0.5))
abline(h=true_ate, lwd=4, col=rgb(0,0.75,0,0.75))
legend('topright', 'True ATE', lwd=4, col=rgb(0,0.75,0,0.75))
plot(0, xlim=range(test_rho), ylim=c(0,max(vars_gp)),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='Estimated Variance',
     main='Robust Variance Estimates')
for(i in 1:n_sims){
    lines(test_rho, vars_gp[i,], lwd=1.5, col=rgb(0,0,0,0.15))
}
lines(test_rho, apply(vars_gp, 2, mean), lwd=4, col=rgb(0.75,0,0,0.75))
lines(test_rho, apply(mean_gp, 2, var), lwd=4, col=rgb(0,0,0.75,0.75))
lines(test_rho, apply(mean_gp, 2, var) /  apply(vars_gp, 2, mean),
      lwd=4, col=rgb(0.75,0,0.75,0.75))

legend('topright', 'Emp. Var.', lwd=4, col=rgb(0,0,0.75,0.75))
plot(0, xlim=range(test_rho), ylim=c(0,1),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='Coverage Proportion',
     main="Coverage")
lines(test_rho, covers, lwd=3, col=rgb(0,0,0,0.5))
abline(h=0.95, lty=3)
abline(h=1, lty=3)

plot(test_rho, apply(mean_gp, 2, var) /  apply(vars_gp, 2, mean),
     type='l', lwd=4, col=rgb(0.75,0,0,0.5),
     xlab='Inverse-Length Scale',
     ylab='Ratio',
     main='Ratio of Emp Var to Mean of Robust Var. Est')
abline(h=1, lwd=4, col=rgb(0,0,0,0.5))

plot(test_rho, apply(ess_gp, 2, mean),
     type='l', lwd=4, col=rgb(0.75,0,0,0.5),
     xlab='Inverse-Length Scale',
     ylab='Ratio',
     main='Ratio of Emp Var to Mean of Robust Var. Est')
abline(h=1, lwd=4, col=rgb(0,0,0,0.5))




# dev.off()



# GLM RESULTS ------------------------------------------------------------------


lowers <- mean_glm - 1.96 * sqrt(vars_glm)
uppers <- mean_glm + 1.96 * sqrt(vars_glm)
covers <- apply((lowers < true_ate) & (uppers > true_ate),2, mean)
plot(test_rho, covers, col=rgb(0,0.75,0.75,0.5), lwd=5)
abline(h=0.95)


plot(0, xlim=range(test_rho), ylim=range(abs(bals_glm)),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='Mean Balance - Abs(Std Diff.)',
     main='Covariate Balance')
for(i in 1:n_sims){
    lines(test_rho, abs(bals_glm[i,]), lwd=1.5, col=rgb(0,0,0,0.15))
}
lines(test_rho, apply(bals_glm, 2, mean),
      lwd=4, col=rgb(0.75,0,0,0.5))
plot(0, xlim=range(test_rho), ylim=range(mean_glm),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='ATE Estimate',
     main='Treatment Effect Estimates')
for(i in 1:n_sims){
    lines(test_rho, mean_glm[i,], lwd=1.5, col=rgb(0,0,0,0.15))
}
lines(test_rho, apply(mean_glm, 2, mean), lwd=4, col=rgb(0.75,0,0,0.5))
abline(h=3, lwd=4, col=rgb(0,0.75,0,0.75))
legend('topright', 'True ATE', lwd=4, col=rgb(0,0.75,0,0.75))
plot(0, xlim=range(test_rho), ylim=c(0,max(vars_glm)),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='Estimated Variance',
     main='Robust Variance Estimates')
for(i in 1:n_sims){
    lines(test_rho, vars_glm[i,], lwd=1.5, col=rgb(0,0,0,0.15))
}
lines(test_rho, apply(vars_glm, 2, mean), lwd=4, col=rgb(0.75,0,0,0.75))
lines(test_rho, apply(mean_glm, 2, var), lwd=4, col=rgb(0,0,0.75,0.75))
lines(test_rho, apply(mean_glm, 2, var) /  apply(vars_glm, 2, mean),
      lwd=4, col=rgb(0.75,0,0.75,0.75))

legend('topright', 'Emp. Var.', lwd=4, col=rgb(0,0,0.75,0.75))
plot(0, xlim=range(test_rho), ylim=c(0,1),
     pch=19, col=rgb(0,0,0,0),
     xlab='Inverse-Length Scale', ylab='Coverage Proportion',
     main="Coverage")
lines(test_rho, covers, lwd=3, col=rgb(0,0,0,0.5))
abline(h=0.95, lty=3)
abline(h=1, lty=3)

plot(test_rho, apply(mean_glm, 2, var) /  apply(vars_glm, 2, mean),
     type='l', lwd=4, col=rgb(0.75,0,0,0.5),
     xlab='Inverse-Length Scale',
     ylab='Ratio',
     main='Ratio of Emp Var to Mean of Robust Var. Est')
abline(h=1, lwd=4, col=rgb(0,0,0,0.5))



lowers <- mean_glm - 1.96 * sqrt(vars_glm)
uppers <- mean_glm + 1.96 * sqrt(vars_glm)
covers <- apply((lowers < true_ate) & (uppers > true_ate),2, mean)
plot(test_rho, covers, ylim=c(0,1), col=rgb(0,0,0.75,0.5), lwd=5, type='l')

lowers1 <- mean_glm - 1.96 * sqrt(ess_glm * vars_glm / n_obs)
uppers1 <- mean_glm + 1.96 * sqrt(ess_glm * vars_glm / n_obs)
covers1 <- apply((lowers1 < true_ate) & (uppers1 > true_ate),2, mean)
lines(test_rho, covers1, ylim=c(0,1), col=rgb(0.75,0,0.75,0.5), lwd=5, type='l')

emp_var <- matrix(apply(mean_glm, 2, var), nrow=n_sims, ncol = length(test_rho), byrow = T)

lowers2 <- mean_glm - 1.96 * sqrt(emp_var)
uppers2 <- mean_glm + 1.96 * sqrt(emp_var)
covers2 <- apply((lowers2 < true_ate) & (uppers2 > true_ate),2, mean)
lines(test_rho, covers2, ylim=c(0,1), col=rgb(0,0.75,0,0.5), lwd=5, type='l')
abline(h=0.95)




# plot(apply(ess2_gp,2,mean)/n_obs, apply(mean_gp, 2, var) /  apply(vars_gp, 2, mean))
# abline(h=1)
# #
# #
# # dset <- data.frame("YO" = yo, "TA" = t_test, "WTS" = wts_old)
# # design.ps <- survey::svydesign(ids=~1, weights=~WTS, data=dset)
# # glm1 <- survey::svyglm(YO ~ TA, design=design.ps)
# # summary(glm1)$coef[2,1:2]
