make_binary <- function(class_labels){
    classes <- as.vector(unique(class_labels))
    classes <- classes[order(classes)]
    n_unique <- length(classes)
    out_mat <- matrix(NA, nrow=length(class_labels), ncol=n_unique)
    for(c in 1:n_unique){
        out_mat[,c] <- ifelse(class_labels == classes[c], 1, 0)
    }
    as.vector(out_mat)
}

mcla_optimize <- function(X,
                          y,
                          cov_function,
                          n_classes,
                          init_theta,
                          return_theta = F,
                          verbose = F){

    objective_function <- function(theta){
        cov_matrix <- cov_function(as.matrix(X), rep(theta, n_classes))
        ps_res <- gp_mcla(cov_matrix, y, n_classes, verbose = F, tol = 1e-6)
        return(-ps_res$logmarglik)
    }

    start_time <- Sys.time()
    if(verbose){
        message(paste('Starting Optimization  @  ', start_time))
    }
    opt_theta <- minqa::bobyqa(par = init_theta,
                               fn = objective_function,
                               lower = rep(0, length(init_theta)),
                               control=list('maxfun'=200))
    print(opt_theta)
    end_time <- Sys.time()
    if(verbose){
        message(paste('Finished Optimization  @  ', end_time))
        message(paste('Time Difference          :', round(difftime(end_time, start_time, units='secs'), 4)))
        message(paste('Optimal Covariate Balance:', opt_theta$fval))
    }

    opt_matrix <- cov_function(as.matrix(X), rep(opt_theta$par, n_classes))
    opt_ps <- gp_mcla(opt_matrix, y, n_classes, verbose = F, tol = 1e-6)
    opt_ps$ComputationTime <- difftime(end_time, start_time, units='secs')
    opt_ps$thetas <- opt_theta$par

    return(opt_ps)

}
