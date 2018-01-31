set.seed(1172018)
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

lm_ps <- function(Y, X, wts, true_val = NULL){
    W <- diag(wts)
    invXtWX <- solve(t(X) %*% W %*% X)
    betas <- invXtWX %*% t(X) %*% W %*% Y

    Yhat <- X %*% betas
    resids <- Y - Yhat
    sighat <- as.double(sum(resids^2) / (length(Y) - 1))

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

lm_ps_wtdsighat <- function(Y, X, wts, true_val = NULL){
    W <- diag(wts)
    invXtWX <- solve(t(X) %*% W %*% X)
    betas <- invXtWX %*% t(X) %*% W %*% Y

    Yhat <- X %*% betas
    resids <- Y - Yhat
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


n_sims <- 1000
true_ate <- 3
mean_true <- matrix(NA, nrow=n_sims, ncol=3)
mean_glmp <- matrix(NA, nrow=n_sims, ncol=3)
vars_true <- matrix(NA, nrow=n_sims, ncol=3)
vars_glmp <- matrix(NA, nrow=n_sims, ncol=3)
covs_true <- matrix(NA, nrow=n_sims, ncol=3)
covs_glmp <- matrix(NA, nrow=n_sims, ncol=3)
ess_true <- matrix(NA, nrow=n_sims, ncol=2)
ess_glmp <- matrix(NA, nrow=n_sims, ncol=2)

true_att <- 5.8125

for(i in 1:n_sims){
    n_obs <- 5000
    Xt <- rnorm(n_obs/5, mean=2.5, sd=0.25)
    Xc <- rnorm(n_obs - n_obs/5, mean=0, sd=2)
    X <- c(Xt, Xc)
    TA <- c(rep(1, n_obs/5), rep(0, n_obs - n_obs/5))
    design_mat <- cbind(1, TA)

    po_stdev <- 0.25
    YT <- X^2 + 2 + rnorm(n_obs, 0, po_stdev)
    YC <- X + rnorm(n_obs, 0, po_stdev)
    YO <- TA * YT + (1-TA)*YC

    est_p <- glm(TA ~ X + I(X^2), family='binomial')
    wts2 <- make_wts2(TA, fitted(est_p))

    ess_glmp[i, ] <- c(sum(wts2)^2 / sum(wts2^2),
                       sum(wts2[TA==1])^2 / sum(wts2[TA==1]^2) +
                           sum(wts2[TA==0])^2 / sum(wts2[TA==0]^2))

    fit_glmp1 <- lm_ps(YO, design_mat, wts=wts2, true_att)
    fit_glmp2 <- lm_ps_wtdsighat(YO, design_mat, wts=wts2, true_att)
    fit_glmp3 <- lm_ps_robust(YO, design_mat, wts=wts2, true_att)

    mean_glmp[i,] <- c(fit_glmp1$ests$coef[2],
                       fit_glmp2$ests$coef[2],
                       fit_glmp3$ests$coef[2])
    vars_glmp[i,] <- c(fit_glmp1$ests$stderrs[2]^2,
                       fit_glmp2$ests$stderrs[2]^2,
                       fit_glmp3$ests$stderrs[2]^2)
    covs_glmp[i,] <- c(fit_glmp1$covers,
                       fit_glmp2$covers,
                       fit_glmp3$covers)
}
apply(covs_true,2,mean)
apply(covs_glmp,2,mean)

lowers1 <- mean_glmp - 1.96 * sqrt(vars_glmp)
uppers1 <- mean_glmp + 1.96 * sqrt(vars_glmp)
covers1 <- apply((lowers1 < true_att) & (uppers1 > true_att), 2, mean)

lowers2 <- mean_glmp - 1.96 * sqrt(ess_glmp[,2] * vars_glmp / n_obs)
uppers2 <- mean_glmp + 1.96 * sqrt(ess_glmp[,2] * vars_glmp / n_obs)
covers2 <- apply((lowers2 < true_att) & (uppers2 > true_att), 2, mean)

print(rbind(covers1, covers2))
