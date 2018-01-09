gpr <- function(X_train,
                y_target,
                cov_func,
                hyperparams,
                noise,
                X_test){
    n_train <- nrow(X_train)
    K_big <- cov_func(rbind(X_train,X_test), hyperparams, noise=0)
    outro = c_gpr(K_big[1:n_train, 1:n_train],
          y_target,
          K_big[1:n_train, (n_train+1):(n_train+nrow(X_test))],
          K_big[(n_train+1):(n_train+nrow(X_test)), (n_train+1):(n_train+nrow(X_test))],
          noise = noise)
    return(outro)
}
