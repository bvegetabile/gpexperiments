
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



library('gpexperiments')
X <- as.matrix(seq(-3, 3, length.out = 100))

cov_sqexp <- gpexperiments::sqexp_common(X, 1, 1, noise = 1e-6)
cov_poly1 <- polykernel(X, 0.5, 1)
cov_poly2 <- polykernel(X, 0.5, 2)

cov_added <- cov_sqexp + cov_poly1

par(mfrow=c(2,3))
image(X, X, cov_sqexp)
image(X, X, cov_poly1)
image(X, X, cov_poly2)
image(X, X, cov_sqexp + cov_poly1)
image(X, X, cov_sqexp + cov_poly2)
image(X, X, cov_sqexp + cov_poly2 + cov_poly1)

n_samples <- 10
draws1 <- mvtnorm::rmvnorm(n_samples, sigma=cov_sqexp)
draws2 <- mvtnorm::rmvnorm(n_samples, sigma=cov_poly1)
draws3 <- mvtnorm::rmvnorm(n_samples, sigma=cov_added)

alldraws <- rbind(draws1, draws2, draws3)

plot(0, xlim=c(-3,3), ylim=range(alldraws),
     xlab='', ylab='', pch=19, col=rgb(0,0,0,0))
for(i in 1:n_samples){
    lines(X, draws1[i,], lwd=3, col=rgb(0,0,0,0.25))
}
plot(0, xlim=c(-3,3), ylim=range(alldraws),
     xlab='', ylab='', pch=19, col=rgb(0,0,0,0))
for(i in 1:n_samples){
    lines(X, draws2[i,], lwd=3, col=rgb(0,0,0,0.25))
}
plot(0, xlim=c(-3,3), ylim=range(alldraws),
     xlab='', ylab='', pch=19, col=rgb(0,0,0,0))
for(i in 1:n_samples){
    lines(X, draws3[i,], lwd=3, col=rgb(0,0,0,0.25))
}


n_test <- 500
X_test <- as.matrix(rnorm(n_test))
p_test <- pnorm(-0.5 * X_test)
t_test <- rbinom(n_test, 1, p_test)
plot(X_test, p_test, ylim=c(0,1))

sqexp_cov <- sqexp_common(X_test, 1, noise=1e-6)
poly1_cov <- polykernel(X_test, 0.75, pwr = 1)
poly2_cov <- polykernel(X_test, 0.75, pwr = 1)

fit1 <- gpbalancer::gpbal_fixed(t_test, sqexp_cov)
fit2 <- gpbalancer::gpbal_fixed(t_test, poly1_cov)
fit3 <- gpbalancer::gpbal_fixed(t_test, poly2_cov)
fit4 <- gpbalancer::gpbal_fixed(t_test, sqexp_cov + poly1_cov)
fit5 <- gpbalancer::gpbal_fixed(t_test, sqexp_cov + poly2_cov)
fit6 <- gpbalancer::gpbal_fixed(t_test, sqexp_cov + poly1_cov + poly2_cov)

par(mfrow=c(2,3))
plot(p_test, fit1$ps)
plot(p_test, fit2$ps)
plot(p_test, fit3$ps)
plot(p_test, fit4$ps)
plot(p_test, fit5$ps)
plot(p_test, fit6$ps)


n_test <- 1500
X_test <- as.matrix(rnorm(n_test))
p_test <- 0.8 * pnorm(-0.5 * X_test^3 + 0.25 * X_test) + 0.1
p_test <- 0.8 * pnorm(0.75 * X_test^2) - 0.25 + 0.1 *X_test
t_test <- rbinom(n_test, 1, p_test)
plot(X_test, p_test, ylim=c(0,1))

sqexp_cov <- sqexp(X_test, c(3), noise=1e-6)
poly1_cov <- polykernel(X_test, 0.5, pwr = 1)
poly2_cov <- polykernel(X_test, 0.5, pwr = 2)

fit1 <- gpbalancer::gpbal_fixed(t_test, sqexp_cov)
fit2 <- gpbalancer::gpbal_fixed(t_test, poly1_cov)
fit3 <- gpbalancer::gpbal_fixed(t_test, poly2_cov)
fit4 <- gpbalancer::gpbal_fixed(t_test, sqexp_cov + poly1_cov)
fit5 <- gpbalancer::gpbal_fixed(t_test, sqexp_cov + poly2_cov)
fit6 <- gpbalancer::gpbal_fixed(t_test, sqexp_cov + poly1_cov + poly2_cov)

par(mfrow=c(2,3))
plot(X_test, fit1$ps, ylim=c(0,1))
lines(X_test[order(X_test)], p_test[order(X_test)],
      lwd=3, col=rgb(0,0,0,0.5))
plot(X_test, fit2$ps, ylim=c(0,1))
lines(X_test[order(X_test)], p_test[order(X_test)],
      lwd=3, col=rgb(0,0,0,0.5))
plot(X_test, fit3$ps, ylim=c(0,1))
lines(X_test[order(X_test)], p_test[order(X_test)],
      lwd=3, col=rgb(0,0,0,0.5))
plot(X_test, fit4$ps, ylim=c(0,1))
lines(X_test[order(X_test)], p_test[order(X_test)],
      lwd=3, col=rgb(0,0,0,0.5))
plot(X_test, fit5$ps, ylim=c(0,1))
lines(X_test[order(X_test)], p_test[order(X_test)],
      lwd=3, col=rgb(0,0,0,0.5))
plot(X_test, fit6$ps, ylim=c(0,1))
lines(X_test[order(X_test)], p_test[order(X_test)],
      lwd=3, col=rgb(0,0,0,0.5))

# plot(fit1$ps, fit6$ps)
# abline(0,1)

testing <- gpexperiments::gpbal_test(X_test, t_test,
                                     sqexp_poly, c(1,1,1),
                                     return_theta = T, verbose = T)
retest <- gpexperiments::gpbal(X_test, t_test,
                               gpbalancer::par_sqexp, c(1),
                               return_theta=T, verbose=T)

par(mfrow=c(1,2))
plot(X_test, testing$ps, ylim=c(0,1))
lines(X_test[order(X_test)], p_test[order(X_test)],
      lwd=3, col=rgb(0,0,0,0.5))
plot(X_test, retest$ps, ylim=c(0,1))
lines(X_test[order(X_test)], p_test[order(X_test)],
      lwd=3, col=rgb(0,0,0,0.5))


make_wts <- function(ta, ps){
    wts <- data.frame(t=(ta/ps) / sum(ta/ps),
                      c=((1-ta)/(1-ps)) / sum((1-ta)/(1-ps)))
    wts <- ifelse(ta==1, wts$t, wts$c)
    return(wts)
}

n_sims <- 100
outres <- matrix(NA, nrow=n_sims, ncol = 6)
outre2 <- matrix(NA, nrow=n_sims, ncol = 6)
for(s in 1:n_sims){
    message(paste(s,':',sep=""), appendLF = F)
    n_test <- 250
    X_test <- as.matrix(rnorm(n_test))
    p_test <- 0.7 * pnorm(0.5 * X_test^3 + 0.25 * X_test) + 0.15
    t_test <- rbinom(n_test, 1, p_test)

    design_x <- cbind(1, t_test)

    yt <- X_test^2 + 2 + rnorm(n_test, sd = 0.25)
    yc <- X_test + rnorm(n_test, sd = 0.25)
    yo <- t_test * yt + (1 - t_test) * yc

    fit_new <- gpexperiments::gpbal_test(X_test, t_test,
                                         # gpbalancer::par_sqexp, c(1,1),
                                         sqexp_poly, c(1,1,1,1),
                                         return_theta = T, verbose = F)
    fit_old <- gpexperiments::gpbal(X_test, t_test,
                                   gpbalancer::par_sqexp, c(1),
                                   return_theta=T, verbose=F)

    wts_true <- make_wts(t_test, p_test)
    wts_new <- make_wts(t_test, fit_new$ps)
    wts_old <- make_wts(t_test, fit_old$ps)

    ate_true <- lm_ps_robust(yo, design_x, wts = wts_true, 3)
    ate_new <- lm_ps_robust(yo, design_x, wts = wts_new, 3)
    ate_old <- lm_ps_robust(yo, design_x, wts = wts_old, 3)

    dset <- data.frame("YO" = yo, "TA" = t_test, "WTS" = wts_old)
    design.ps <- survey::svydesign(ids=~1, weights=~WTS, data=dset)
    glm1 <- survey::svyglm(YO ~ TA, design=design.ps)
    summary(glm1)$coef[2,1:2]

    outres[s, ] <- c(ate_true$ests$coef[2],ate_true$ests$stderrs[2],
                     ate_new$ests$coef[2],ate_new$ests$stderrs[2],
                     ate_old$ests$coef[2],ate_old$ests$stderrs[2])


}
message(paste('Fin',sep=""), appendLF = T)


mean_res <- outres[,c(1,3,5)]
vars_res <- outres[,c(2,4,6)]
apply(mean_res, 2, var) / apply(vars_res, 2, mean)

