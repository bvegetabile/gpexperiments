
gpR <- function(X_train,
                targets,
                cov_func,
                thetas,
                noise_level = 1e-6,
                X_test = NULL){
    n_train <- nrow(X_train)
    if(is.null(X_test)){
        X <- X_train
        K <- par_sqexp(X, thetas)
        L <- t(chol(K + diag(noise_level, nrow(K), ncol(K))))
        alpha <- solve(t(L), solve(L, targets))
        return(K%*%alpha)
    } else {
        X <- rbind(X_train, X_test)
        K_big <- par_sqexp(X, thetas)
        K <- K_big[1:nrow(X_train), 1:nrow(X_train)]
        K_star <- K_big[1:nrow(X_train), (nrow(X_train)+1):(nrow(X_train)+nrow(X_test))]
        L <- t(chol(K + diag(noise_level, nrow(K), ncol(K))))
        alpha <- solve(t(L), solve(L, targets))

        train_mean <- t(K)%*%alpha
        test_mean <- t(K_star)%*%alpha

        v <- solve(L, K_star)

        Var_f <- K_big[(nrow(X_train)+1):(nrow(X_train)+nrow(X_test)),
                       (nrow(X_train)+1):(nrow(X_train)+nrow(X_test))] - t(v)%*%v

        log_marg <- -0.5 * t(targets) %*% alpha - sum(log(diag(L))) - 0.5 * n_train * log(2 * pi)

        res <- list('test mean' = test_mean,
                    'log-likelihood' = log_marg)
        # print(log_marg)
        return(res)
    }

}
