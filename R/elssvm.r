
#%%
library("methods")
library("stats")
library("quadprog")
library("lpSolve")
library("parallel")
library("mlbench")
library("HTLR")
#%%


cstep.elssvm = function(x = NULL, y = NULL, valid_x = NULL, valid_y = NULL, nfolds = 5,
                        lambda_seq = 2^seq(-10, 10, length.out = 100), 
                        kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian", "gaussian2"), kparam = 1,
                        theta = NULL, scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
{
  call = match.call()
  kernel = match.arg(kernel)
  criterion = match.arg(criterion)
  
  if((criterion != "0-1") && (criterion != "hinge"))
  {
    cat("ERROR: Only 0-1 and hinge can be used as criterion!", "\n")
    return(NULL)
  }
  
  out = list()
  p = NCOL(x)
  
  lambda_seq = as.numeric(lambda_seq)
  kparam = as.numeric(kparam)
  
  
  lambda_seq = sort(lambda_seq, decreasing = FALSE)
  
  if (!is.null(valid_x) & !is.null(valid_y)) {
    
    anova_K = make_anovaKernel(x, x, kernel, kparam)
    if (is.null(theta)) {
      theta = rep(1, anova_K$numK)
    }
    K = combine_kernel(anova_K, theta)
    
    anova_K_valid = make_anovaKernel(valid_x, x, kernel, kparam)
    K_valid = combine_kernel(anova_K_valid, theta)
    
    fold_err = mclapply(1:length(lambda_seq),
                        function(j) {
                          error = try({
                            svm_fit = svm_compact(K = K, y = y, lambda = lambda_seq[j], ...)
                          })
                          
                          if (!inherits(error, "try-error")) {
                            pred_val = predict.svm_compact(svm_fit, K_valid)
                            if (criterion == "0-1") {
                              acc = sum(valid_y == pred_val$class) / length(valid_y)
                              err = 1 - acc
                            } else {
                              
                            }
                          } else {
                            svm_fit = NULL
                            err = Inf
                          }
                          return(list(error = err, fit_model = svm_fit))
                        }, mc.cores = nCores)
    
    valid_err = sapply(fold_err, "[[", "error")
    model_list = lapply(fold_err, "[[", "fit_model")
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_param = c(lambda = lambda_seq[opt_ind])
    opt_valid_err = min(valid_err)
    
  } else {
    
    ran = data_split(y, nfolds)
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(lambda_seq))
    
    for (i_cv in 1:nfolds) {
      
      omit = ran == i_cv
      train_x = x[!omit, ]
      train_y = y[!omit]
      valid_x = x[omit, ]
      valid_y = y[omit]
      
      subanova_K = make_anovaKernel(train_x, train_x, kernel, kparam)
      if (is.null(theta)) {
        theta = rep(1, subanova_K$numK)
      }
      subK = combine_kernel(subanova_K, theta)
      
      subanova_K_valid = make_anovaKernel(valid_x, train_x, kernel, kparam)
      subK_valid = combine_kernel(subanova_K_valid, theta)
      
      fold_err = mclapply(1:length(lambda_seq),
                          function(j) {
                            error = try({
                              svm_fit = svm_compact(K = subK, y = train_y, lambda = lambda_seq[j], ...)
                            })
                            
                            if (!inherits(error, "try-error")) {
                              pred_val = predict.svm_compact(svm_fit, subK_valid)
                              if (criterion == "0-1") {
                                acc = sum(valid_y == pred_val$class) / length(valid_y)
                                err = 1 - acc
                              } else {
                                
                              }
                            } else {
                              svm_fit = NULL
                              err = Inf
                            }
                            return(list(error = err, fit_model = svm_fit))
                          }, mc.cores = nCores)
      
      valid_err = sapply(fold_err, "[[", "error")
      valid_err_mat[i_cv, ] = valid_err
    }
    valid_err = colMeans(valid_err_mat)
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_param = c(lambda = lambda_seq[opt_ind])
    opt_valid_err = min(valid_err)
    valid_x = NULL
    valid_y = NULL
  }
  
  out$opt_param = opt_param
  out$opt_valid_err = opt_valid_err
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  out$x = x
  out$y = y
  out$theta = theta
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$kernel = kernel
  out$kparam = kparam
  out$scale = scale
  out$criterion = criterion
  out$nfolds = nfolds
  if (optModel) {
    
    if (!is.null(valid_x) & !is.null(valid_y)) {
      out$opt_model = model_list[[opt_ind]]
    } else {
      anova_K = make_anovaKernel(x, x, kernel, kparam)
      if (is.null(theta)) {
        theta = rep(1, anova_K$numK)
      }
      K = combine_kernel(anova_K, theta)
      opt_model = svm_compact(K = K, y = y, lambda = opt_param["lambda"], ...)
      out$opt_model = opt_model
    }
  }
  out$call = call
  class(out) = "ssvm"
  return(out)
}



thetastep.elssvm = function(object, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, isCombined = TRUE, cv_type = "original", nCores = 1, gamma = gamma,...)
{
  call = match.call()
  out = list()
  lambda_theta_seq = sort(as.numeric(lambda_theta_seq), decreasing = FALSE)
  lambda = object$opt_param["lambda"]
  criterion = object$criterion
  if((criterion != "0-1") && (criterion != "hinge"))
  {
    cat("Only 0-1 and hinge can be used as criterion!", "\n")
    return(NULL)
  }
  kernel = object$kernel
  kparam = object$kparam
  x = object$x
  y = object$y
  # theta = object$theta
  valid_x = object$valid_x
  valid_y = object$valid_y
  
  nfolds = object$nfolds
  
  anova_K = make_anovaKernel(x, x, kernel, kparam)
  
  if (is.null(object$opt_model)) {
    K = combine_kernel(anova_K, object$theta)
    opt_model = svm_compact(K = K, y = y, lambda = lambda, ...)
  } else {
    opt_model = object$opt_model
  }
  
  if (!is.null(valid_x) & !is.null(valid_y)) {
    
  } else {
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(lambda_theta_seq))
    ran = data_split(y, nfolds)
    
    for (i_cv in 1:nfolds) {
      omit = ran == i_cv
      train_x = x[!omit, ]
      train_y = y[!omit]
      valid_x = x[omit, ]
      valid_y = y[omit]
      
      subanova_K = make_anovaKernel(train_x, train_x, kernel, kparam)
      subK = combine_kernel(subanova_K, object$theta)
      subanova_K_valid = make_anovaKernel(valid_x, train_x, kernel, kparam)
      
      init_model = svm_compact(K = subK, y = train_y, lambda = lambda, ...)
      alpha = init_model$alpha
      bias = init_model$bias
      
      fold_err = mclapply(1:length(lambda_theta_seq),
                          function(j) {
                            error = try({
                              theta = findtheta.elssvm(y = train_y, anova_kernel = subanova_K,
                                                    alpha = alpha, bias = bias, lambda = lambda, lambda_theta = lambda_theta_seq[j], gamma = gamma)
                              if (isCombined) {
                                subK = combine_kernel(subanova_K, theta)
                                init_model = svm_compact(K = subK, y = train_y, lambda = lambda, ...)
                              }
                            })
                            
                            if (!inherits(error, "try-error")) {
                              subK_valid = combine_kernel(subanova_K_valid, theta)
                              pred_val = predict.svm_compact(init_model, newK = subK_valid)
                              
                              if (criterion == "0-1") {
                                acc = sum(valid_y == pred_val$class) / length(valid_y)
                                err = 1 - acc
                              } else {
                                
                              }
                            } else {
                              err = Inf
                              theta = rep(0, anova_K$numK)
                            }
                            return(list(error = err, theta = theta))
                          }, mc.cores = nCores)
    valid_err_mat[i_cv, ] = sapply(fold_err, "[[", "error") 
    }
    valid_err = colMeans(valid_err_mat)
    
    if (cv_type == "original") {
      opt_ind = max(which(valid_err == min(valid_err)))
    } else {
      cv_se = (apply(valid_err_mat, 2, sd) / sqrt(nfolds))
      opt_ind = max(which(valid_err == min(valid_err)))
      opt_ind = max(which(valid_err <= (min(valid_err) + cv_se[opt_ind])))
    }
    opt_lambda_theta = lambda_theta_seq[opt_ind]
    opt_valid_err = min(valid_err)
    
    theta_seq_list = mclapply(1:length(lambda_theta_seq),
                              function(j) {
                                error = try({
                                  theta = findtheta.elssvm(y = y, anova_kernel = anova_K, 
                                                        alpha = opt_model$alpha, bias = opt_model$bias, 
                                                        lambda = lambda, lambda_theta = lambda_theta_seq[j], gamma = gamma)
                                })
                                if (inherits(error, "try-error")) {
                                  theta = rep(0, anova_K$numK)
                                }
                                return(theta)
                              }, mc.cores = nCores)
    theta_seq = do.call(cbind, theta_seq_list)
    opt_theta = theta_seq[, opt_ind]
  }
  
  out$opt_lambda_theta = opt_lambda_theta
  out$opt_ind = opt_ind
  out$opt_theta = opt_theta
  out$theta_seq = theta_seq
  out$opt_valid_err = opt_valid_err
  out$vliad_err = valid_err
  return(out)
}


findtheta.elssvm = function(y, anova_kernel, alpha, bias, lambda, lambda_theta, gamma)
{
    if (anova_kernel$numK == 1)
    {
    cat("Only one kernel", "\n")
    return(c(1))
    }
    
    y_temp = factor(y)
    classname = levels(y_temp)
    n_class = length(classname)
    
    y_int = integer(length(y))
    for (j in 1:n_class) {y_int[which(y_temp %in% classname[j])] = j}
    if (is(y, "numeric")) {classname = as.numeric(classname)}
    
    y_int = ifelse(y_int == 1, 1, -1)
    n = length(y_int)

    cvec = alpha * y_int

    dvec = NULL
    A_mat = NULL

    for (i in 1:anova_kernel$numK) {
        temp_d = (1 / 2) * lambda * t(cvec) %*% anova_kernel$K[[i]] %*% cvec + lambda_theta * gamma
        temp_A = y_int * as.vector(anova_kernel$K[[i]] %*% cvec)
    
    if (i == 1) {
        dvec = temp_d
        A_mat = temp_A
    } else {
        dvec = c(dvec, temp_d)
        A_mat = rbind(A_mat, temp_A)
    }
    }


    dvec = -c(dvec, rep(1 / n, n))
    A_mat = cbind(A_mat, diag(1, anova_kernel$numK))
    # A_mat = rbind(A_mat, diag(1, ncol(A_mat)))
    A_theta = cbind(diag(1, n, n), matrix(0, n, anova_kernel$numK))
    A_mat = rbind(A_mat, A_theta)

    bvec = c(1 - bias * y_int, rep(0, anova_kernel$numK)) # theta >= 0

    
    Dmat = NULL
    D_1 = cbind(diag((1 - gamma) * lambda_theta, anova_kernel$numK), matrix(0, anova_kernel$numK, n))
    D_2 = cbind(matrix(0, n, anova_kernel$numK), diag(1/n, n, n))
    Dmat = 2 * rbind(D_1, D_2)
    
    # find solution by QP
    QP = solve.QP(Dmat = Dmat, dvec = dvec, Amat = A_mat, bvec = bvec, meq = 0, factorized = FALSE)$solution

    # find the theta vector only from the solution
    theta = QP[1:anova_kernel$numK]
    return(round(theta, 6))
}


elssvm = function(x = NULL, y = NULL, valid_x = NULL, valid_y = NULL, nfolds = 5,
                lambda_seq = 2^seq(-10, 10, length.out = 100), 
                lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, cv_type = "original", isCombined = TRUE,
                kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian", "gaussian2"), kparam = 1,
                scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, verbose = 1, epsilon_D = 1e-4, gamma = gamma, ...)
{
    out = list()
    cat("Fit c-step \n")
    cstep_fit = cstep.elssvm(x = x, y = y, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds, 
                        lambda_seq = lambda_seq, kernel = kernel, kparam = kparam, theta = NULL,
                        scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, epsilon_D = epsilon_D, ...)

    cat("Fit theta-step \n")
    thetastep_fit = thetastep.elssvm(cstep_fit, lambda_theta_seq = lambda_theta_seq, isCombined = isCombined, cv_type = "original",  nCores = nCores, epsilon_D = epsilon_D, gamma = gamma, ...)

    cat("Fit c-step \n")
    opt_cstep_fit = cstep.elssvm(x = x, y = y, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds, 
                        lambda_seq = lambda_seq, kernel = kernel, kparam = kparam, theta = thetastep_fit$opt_theta,
                        scale = scale, criterion = criterion, optModel = TRUE, nCores = nCores, epsilon_D = epsilon_D,...)

    if (verbose == 1) {
        cat("CV-error(cstep):", cstep_fit$opt_valid_err, "\n")
        cat("CV-error(theta-step):", thetastep_fit$opt_valid_err, "\n")
        cat("CV-error(cstep):", opt_cstep_fit$opt_valid_err, "\n")
    }

    out$opt_param = opt_cstep_fit$opt_param
    out$x = opt_cstep_fit$x
    out$y = opt_cstep_fit$y
    out$theta = opt_cstep_fit$theta
    out$kernel = opt_cstep_fit$kernel
    out$kparam = opt_cstep_fit$kparam
    out$scale = opt_cstep_fit$scale
    out$criterion = opt_cstep_fit$criterion
    out$nfolds = opt_cstep_fit$nfolds
    out$opt_model = opt_cstep_fit$opt_model
    return(out)
}



predict.elssvm = function(object, newx = NULL)
{   
    model = object$opt_model
    K_valid = combine_kernel(make_anovaKernel(newx, object$x, object$kernel, object$kparam), object$theta)
    pred_y = predict.svm_compact(model, newK = K_valid)$pred_value
    pred_y[pred_y == 0] = 1e-10
    pred_class = sign(pred_y)
    pred_class = ifelse(pred_class == 1, model$classname[1], model$classname[2])

    if (is(object$y, "factor")) {pred_class = factor(pred_class, model$classname)}
    if (is(object$y, "numeric")) {pred_class = as.numeric(pred_class)}
    
    return(list(class = pred_class, pred_value = pred_y))
}




