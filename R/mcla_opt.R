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
        cov_matrix <- cov_function(as.matrix(X), theta)
        ps_res <- gp_mcla(cov_matrix, y, verbose = F, tol = 1e-2)
        ps_est <- pnorm(ps_res$PosteriorMean)

        if(tolower(wts_vers) =='ate'){
            ps_wts <- ifelse(y==1, 1/ps_est, 1/(1-ps_est))
        } else if(tolower(wts_vers) == 'att') {
            ps_wts <- ifelse(y==1, 1, ps_est/(1-ps_est))
        } else {
            message('invalid weighting scheme')
            return(NULL)
        }

        if(balance_metric == 'mom_sq'){
            cb_bal <- .mom_sq_bal(data.frame(X), 1:ncol(X), y==1, ps_wts)
        } else if(balance_metric == 'mom'){
            cb_bal <- .mom_bal(data.frame(X), 1:ncol(X), y==1, ps_wts)
        } else if(balance_metric == 'ks'){
            cb_bal <- 0
            for(i in 1:ncol(X)){
                cb_bal <- cb_bal + .ks_avg_test(X[,i], y, ps_est, 500)
            }
        } else(
            return(NULL)
        )
        return(cb_bal)
    }

}
