
library('gpexperiments')

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

make_wts <- function(ta, ps){
    wts <- data.frame(t=ta / sum(ta),
                      c=((1-ta) * ps /(1-ps)) / sum((1-ta) * ps/(1-ps)))
    wts <- ifelse(ta==1, wts$t, wts$c)
    return(wts)
}

n_sims <- 1
outres <- matrix(NA, nrow=n_sims, ncol = 6)
balres <- matrix(NA, nrow=n_sims, ncol = 4)
outre2 <- matrix(NA, nrow=n_sims, ncol = 6)
for(s in 1:n_sims){
    message(paste(s,':',sep=""), appendLF = F)
    n_test <- 250
    Xt <- rnorm(n_test/5, mean=2.5, sd=0.25)
    Xc <- rnorm(n_test - n_test/5, mean=0, sd=2)
    X_test <- as.matrix(c(Xt, Xc))
    t_test <- c(rep(1, n_test/5), rep(0, n_test - n_test/5))

    p_test <- fitted(glm(t_test ~ X_test + I(X_test^2), family='binomial'))

    design_x <- cbind(1, t_test)

    yt <- X_test^2 + 2 + rnorm(n_test, sd = 0.25)
    yc <- X_test + rnorm(n_test, sd = 0.25)
    yo <- t_test * yt + (1 - t_test) * yc

    fit_new <- gpexperiments::gpbal_test(X_test, t_test,
                                         sqexp_poly, c(1,1,1),
                                         return_theta = T, verbose = F)
    fit_old <- gpexperiments::gpbal(X_test, t_test,
                                   gpbalancer::par_sqexp, c(1),
                                   return_theta=T, verbose=F)

    wts_true <- make_wts(t_test, p_test)
    wts_new <- make_wts(t_test, fit_new$ps)
    wts_old <- make_wts(t_test, fit_old$ps)

    bal_b4 <- gpbalancer::bal_stats(X_test, t_test, 'continuous')[7]
    bal_tru <- gpbalancer::bal_stats(X_test, t_test, 'continuous', wts = wts_true)[7]
    bal_new <- gpbalancer::bal_stats(X_test, t_test, 'continuous', wts = wts_new)[7]
    bal_old <- gpbalancer::bal_stats(X_test, t_test, 'continuous', wts = wts_old)[7]

    att_true <- lm_ps_robust(yo, design_x, wts = wts_true, 5.8125)
    att_new <- lm_ps_robust(yo, design_x, wts = wts_new, 5.8125)
    att_old <- lm_ps_robust(yo, design_x, wts = wts_old, 5.8125)

    dset <- data.frame("YO" = yo, "TA" = t_test, "WTS" = wts_old)
    design.ps <- survey::svydesign(ids=~1, weights=~WTS, data=dset)
    glm1 <- survey::svyglm(YO ~ TA, design=design.ps)
    summary(glm1)$coef[2,1:2]

    outres[s, ] <- c(att_true$ests$coef[2],att_true$ests$stderrs[2],
                     att_new$ests$coef[2],att_new$ests$stderrs[2],
                     att_old$ests$coef[2],att_old$ests$stderrs[2])
    balres[s, ] <- c(bal_b4, bal_tru, bal_new, bal_old)
}
message(paste('Fin',sep=""), appendLF = T)

pdf('comparing_kernel_estimates.pdf', height=3.5, width=10)
par(mfrow=c(1,2))
plot(X_test[order(X_test)], fit_old$ps[order(X_test)],
     type='l',
     lty=1, lwd=3, col=rgb(0.75,0,0,0.5),
     ylim=c(0,1),
     xlab='Covariate Value',
     ylab='P(t=1|x)',
     main='Squared Exponential Kernel Only')
plot(X_test[order(X_test)], fit_new$ps[order(X_test)],
     type='l',
     lty=1, lwd=3, col=rgb(0.25,0,0.75,0.5),
     ylim=c(0,1),
     xlab='Covariate Value',
     ylab='P(t=1|x)',
     main='Squared Exp + Normalized Poly (p=1) Kernel')
dev.off()



mean_res <- outres[,c(1,3,5)]
vars_res <- outres[,c(2,4,6)]
apply(mean_res, 2, var) / apply(vars_res, 2, mean)

