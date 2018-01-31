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


n_sims <- 100
true_ate <- 3
mean_true <- matrix(NA, nrow=n_sims, ncol=3)
mean_glmp <- matrix(NA, nrow=n_sims, ncol=3)
vars_true <- matrix(NA, nrow=n_sims, ncol=3)
vars_glmp <- matrix(NA, nrow=n_sims, ncol=3)
covs_true <- matrix(NA, nrow=n_sims, ncol=3)
covs_glmp <- matrix(NA, nrow=n_sims, ncol=3)
ess_true <- matrix(NA, nrow=n_sims, ncol=2)
ess_glmp <- matrix(NA, nrow=n_sims, ncol=2)

for(i in 1:n_sims){
    n_obs <- 500
    X <- rnorm(n_obs)
    p <- pnorm(-0.5 * X)
    TA <- rbinom(n_obs, 1, p)
    design_mat <- cbind(1, TA)

    YT <- X^2 + 2 + rnorm(n_obs, 0, 0.25)
    YC <- X + rnorm(n_obs, 0, 2)
    YO <- TA * YT + (1-TA)*YC

    VarYO <- var(YO)

    est_p <- glm(TA ~ X, family='binomial')
    wts1 <- make_wts(TA, p)
    wts2 <- make_wts(TA, fitted(est_p))

    ess_test <- (sum(wts[TA==1])*sum(wts[TA==0]))^(-2) * sum( (TA * sum(wts) - sum(wts[TA==1]))^2)


    ess_true[i, ] <- c(sum(wts1)^2 / sum(wts1^2),
                       sum(wts1[TA==1])^2 / sum(wts1[TA==1]^2) +
                           sum(wts1[TA==0])^2 / sum(wts1[TA==0]^2))
    ess_glmp[i, ] <- c(sum(wts2)^2 / sum(wts2^2),
                       sum(wts2[TA==1])^2 / sum(wts2[TA==1]^2) +
                           sum(wts2[TA==0])^2 / sum(wts2[TA==0]^2))

    fit_true1 <- lm_ps(YO, design_mat, wts=wts1, 3)
    fit_true2 <- lm_ps_wtdsighat(YO, design_mat, wts=wts1, 3)
    fit_true3 <- lm_ps_robust(YO, design_mat, wts=wts1, 3)

    fit_glmp1 <- lm_ps(YO, design_mat, wts=wts2, 3)
    fit_glmp2 <- lm_ps_wtdsighat(YO, design_mat, wts=wts2, 3)
    fit_glmp3 <- lm_ps_robust(YO, design_mat, wts=wts2, 3)

    mean_true[i,] <- c(fit_true1$ests$coef[2],
                       fit_true2$ests$coef[2],
                       fit_true3$ests$coef[2])
    vars_true[i,] <- c(fit_true1$ests$stderrs[2]^2,
                       fit_true2$ests$stderrs[2]^2,
                       fit_true3$ests$stderrs[2]^2)
    covs_true[i,] <- c(fit_true1$covers,
                       fit_true2$covers,
                       fit_true3$covers)

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
covers1 <- apply((lowers1 < 3) & (uppers1 > 3), 2, mean)

lowers2 <- mean_glmp - 1.96 * sqrt(ess_glmp[,2] * vars_glmp / n_obs)
uppers2 <- mean_glmp + 1.96 * sqrt(ess_glmp[,2] * vars_glmp / n_obs)
covers2 <- apply((lowers2 < 3) & (uppers2 > 3), 2, mean)

print(rbind(covers1, covers2))
