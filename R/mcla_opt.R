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


multicovbal <- function(X, TA, wts=NULL, plot=T){
    X <- as.matrix(X)
    lvls <- unique(TA)
    n_lvls <- length(lvls)

    n_covs <- ncol(X)
    barx_cov <- apply(X, 2, mean)
    varx_cov <- apply(X, 2, var)

    before_table <- matrix(nrow = n_covs, ncol = 2*n_lvls + 1)
    colnames(before_table)[1:(2*n_lvls + 1)] <- c('Overall Mean')
    colnames(before_table)[seq(2, 2*n_lvls, 2)] <- paste("TA=", lvls, ":Mean", sep = '')
    colnames(before_table)[seq(3, 2*n_lvls+1, 2)] <- paste("TA=", lvls, ":Delta", sep = '')
    for(i in 1:n_covs){
        before_table[i, 1] <- barx_cov[i]
        for(j in 1:n_lvls){
            barx_cov_lvl <- mean(X[TA==lvls[j], i])
            varx_cov_lvl <- var(X[TA==lvls[j], i])
            before_table[i, 2*j] <- barx_cov_lvl
            before_table[i, 2*j + 1] <- (barx_cov[i] - barx_cov_lvl) / sqrt((varx_cov[i] + varx_cov_lvl)/2)
        }
    }

    if(is.null(wts)){
        return(before_table)
    } else {
        after_table <- matrix(nrow = n_covs, ncol = 2*n_lvls + 1)
        colnames(after_table)[1:(2*n_lvls + 1)] <- c('Overall Mean')
        colnames(after_table)[seq(2, 2*n_lvls, 2)] <- paste("TA=", lvls, ":Mean", sep = '')
        colnames(after_table)[seq(3, 2*n_lvls+1, 2)] <- paste("TA=", lvls, ":Delta", sep = '')
        for(i in 1:n_covs){
            after_table[i, 1] <- barx_cov[i]
            for(j in 1:n_lvls){
                barx_cov_lvl <- sum(X[TA==lvls[j] ,i] * wts[TA==lvls[j]]) / sum(wts[TA==lvls[j]])
                varx_cov_lvl <- sum((X[TA==lvls[j] ,i] - barx_cov_lvl)^2 * wts[TA==lvls[j]]) / sum(wts[TA==lvls[j]])
                after_table[i, 2*j] <- barx_cov_lvl
                after_table[i, 2*j + 1] <- (barx_cov[i] - barx_cov_lvl) / sqrt((varx_cov[i] + varx_cov_lvl)/2)
            }
        }
        return(list('before_balance' = round(before_table,3),
                    'after_balance' = round(after_table,3)))
    }
}



mcgpbal <- function(){

}
