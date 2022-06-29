#%%
library("methods")
library("stats")
library("quadprog")
library("lpSolve")
library("parallel")
library("mlbench")
library("HTLR")
#%%

train_data = bsim_gen(50, 100, seed = 1, type = "neuralnet")
train_x = scale(train_data$x)
train_y = as.numeric(as.character(train_data$y))

test_data = bsim_gen(10000, 100, seed = 1, type = "neuralnet")
test_x = scale(test_data$x)
test_y = as.numeric(as.character(test_data$y))

############################################################################
# cosso - svm

cstep.out = cstep.ssvm(x = train_x, y = train_y, nfolds = 5, lambda_seq = 2^{seq(-10, 10, length.out = 100)}, kernel ="linear", kparam = 1, theta = NULL, scale = FALSE,
                        criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, epsilon_D = 1e-4)

thetastep.out = thetastep.ssvm(object = cstep.out, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, isCombined = TRUE, cv_type = "original", nCores = 1, epsilon_D = 1e-4)


c1step.out = cstep.ssvm(x = train_x, y = train_y, nfolds = 5, lambda_seq = 2^seq(-10, 10, length.out = 100), kernel ="linear", kparam = 1, theta = thetastep.out$opt_theta, scale = FALSE,
                        criterion = c("0-1", "loss"), optModel = TRUE, nCores = 1, epsilon_D = 1e-4)

ssvm_result = ssvm(x = train_x, y = train_y, nfolds = 5,lambda_seq = 2^seq(-10, 10, length.out = 100), 
                lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, cv_type = "original", isCombined = TRUE,
                kernel = "spline", kparam = 1,
                scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, verbose = 1, epsilon_D = 1)

result = predict.ssvm(object = ssvm_result, newx = test_x)
score0 = score(y_test = test_y, pred = result$class)

#############################################################################
# elcosso - svm

elcstep.out = cstep.elssvm(x = train_x, y = train_y, nfolds = 5, lambda_seq = 2^{seq(-10, 10, length.out = 100)}, kernel ="linear", kparam = 1, theta = NULL, scale = FALSE,
                        criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, epsilon_D = 1e-4)

elthetastep.out = thetastep.elssvm(object = elcstep.out, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, isCombined = TRUE, cv_type = "original", nCores = 1, epsilon_D = 1e-4, gamma = 0.5)


elc1step.out = cstep.elssvm(x = train_x, y = train_y, nfolds = 5, lambda_seq = 2^seq(-10, 10, length.out = 100), kernel ="linear", kparam = 1, theta = elthetastep.out$opt_theta, scale = FALSE,
                        criterion = c("0-1", "loss"), optModel = TRUE, nCores = 1, epsilon_D = 1e-1)

elresult = predict.elssvm(object = elc1step.out, newx = test_x)
score_elssvm = score(y_test = test_y, pred = elresult$class)

elssvm_result = elssvm(x = train_x, y = train_y, nfolds = 5,lambda_seq = 2^seq(-10, 10, length.out = 100), 
                lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, cv_type = "original", isCombined = TRUE,
                kernel = "spline", kparam = 1,
                scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, verbose = 1, epsilon_D = 1, gamma = 0.5)

elresult = predict.elssvm(object = elssvm_result, newx = test_x)
score1 = score(y_test = test_y, pred = elresult$class)



###########################################################################
#2022 - 03 - 00 ~ simulation, linear beta() data simulation
###########################################################################

opt_ind = c()

ssvm_f1 = c()
ssvm_er = c()

elssvm_f1 = c()
elssvm_er = c()

grid_f1 = c()
grid_er = c()

gamma = 2^{seq(-2, -0.1e-3, length.out = 5)}


for (i in 1:100) {
cat("simulation_number:", i, "\n")

train_data = bsim_gen(50, 200, seed = i, type = "neuralnet")
train_x = scale(train_data$x)
train_y = as.numeric(as.character(train_data$y))

test_data = bsim_gen(10000, 200, seed = i, type = "neuralnet")
test_x = scale(test_data$x)
test_y = as.numeric(as.character(test_data$y))

ssvm_out = ssvm(x = train_x, y = train_y, nfolds = 5,lambda_seq = 2^seq(-10, 10, length.out = 100), 
                lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, cv_type = "original", isCombined = TRUE,
                kernel = "spline", kparam = 1,
                scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, verbose = 1, epsilon_D = 1)


ssvm_pred = predict.ssvm(object = ssvm_out, newx = test_x)
ssvmresult = score(y_test = test_y, pred = ssvm_pred$class)

ssvm_f1[i] = ssvmresult$f1score
ssvm_er[i] = ssvmresult$error_rate  

if (ssvmresult$Recall == 0){

elssvm_f1[i] = NaN
elssvm_er[i] = NaN

} else {
    for (j in 1:5){
        
            elssvm_out = elssvm(x = train_x, y = train_y, nfolds = 5,lambda_seq = 2^seq(-10, 10, length.out = 100), 
                    lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, cv_type = "original", isCombined = TRUE,
                    kernel = "spline", kparam = 1,
                    scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, verbose = 1, epsilon_D = 1, gamma = gamma[j])

    grid_elssvm_pred = predict.ssvm(object = elssvm_out, newx = test_x)
    grid_elssvmresult = score(y_test = test_y, pred = grid_elssvm_pred$class)
    
    grid_f1[j] = grid_elssvmresult$f1score
    grid_er[j] = grid_elssvmresult$error_rate
    }


opt_ind = which(grid_f1 == max(grid_f1))[1]

elssvm_f1[i] = grid_f1[opt_ind]
elssvm_er[i] = grid_er[opt_ind]
}
}



ind = which(is.na(ssvm_f1) == TRUE | is.na(elssvm_f1) == TRUE)
ssvm_f1 = ssvm_f1[-c(ind)]
elssvm_f1 = elssvm_f1[-c(ind)]
ssvm_sd = ssvm_sd[-c(ind)]
elssvm_sd = elssvm_sd[-c(ind)]


mean(na.omit(ssvm_f1))
mean(na.omit(elssvm_f1))

mean(na.omit(ssvm_er))
mean(na.omit(elssvm_er))

sd(na.omit(ssvm_er))
sd(na.omit(elssvm_er))



r = 2

    rho = matrix(0, r, r)
    
    diag(rho) = 1

    for (i in 1:2){
        for (j in 1:2){
            if (i == j){
                rho[i,j] = 1
            }
            else{
                rho[i,j] = 0.5
            }
    }
    }

    mutemp = rep(0, 2)
    
    X_tmp = mvrnorm(n = 10, mu = mutemp, Sigma = rho)
    X_tmp = cbind(X_tmp, X_tmp[, 1] * X_tmp[, 2])

