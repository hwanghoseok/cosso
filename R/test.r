
############################################################################


train_data = bsim_gen(50, 10, seed = 1, type = "neralnet")
train_x = scale(train_data$x)
train_y = as.numeric(as.character(train_data$y))

test_data = bsim_gen(10000, 10, seed = 1, type = "neralnet")
test_x = scale(test_data$x)
test_y = as.numeric(as.character(test_data$y))

############################################################################
# cosso - svm

cstep.out = cstep.ssvm(x = train_x, y = train_y, nfolds = 5, lambda_seq = 2^{seq(-10, 10, length.out = 100)}, kernel ="gaussian", kparam = 1, theta = NULL, scale = FALSE,
                        criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, epsilon_D = 1e-4)

thetastep.out = thetastep.ssvm(object = cstep.out, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, isCombined = TRUE, cv_type = "original", nCores = 1, epsilon_D = 1e-4)


c1step.out = cstep.ssvm(x = train_x, y = train_y, nfolds = 5, lambda_seq = 2^seq(-10, 10, length.out = 100), kernel ="gaussian", kparam = 1, theta = thetastep.out$opt_theta, scale = FALSE,
                        criterion = c("0-1", "loss"), optModel = TRUE, nCores = 1, epsilon_D = 1e-4)

ssvm_result = ssvm(x = train_x, y = train_y, nfolds = 5,lambda_seq = 2^seq(-10, 10, length.out = 100), 
                lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, cv_type = "original", isCombined = TRUE,
                kernel = "gaussian", kparam = 1,
                scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, verbose = 1, epsilon_D = 1e-6)

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
                kernel = "gaussian", kparam = 1,
                scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, verbose = 1, epsilon_D = 1e-4, gamma = 0.9)

elresult = predict.ssvm(object = ssvm_result, newx = test_x)
score1 = score(y_test = test_y, pred = elresult$class)

#################################################################################
# 변수선택, 성능비교

true_ind = which(train_data$true == 1)

ssvm_ind = which(ssvm_result$theta != 0)
ssvm_select = intersect(true_ind, ssvm_ind) # true 변수중에 선택된거

ssvm_recall = length(ssvm_select)/length(true_ind) # 실제 true중 선택된 변수 비율
ssvm_precision = length(ssvm_select) / length(ssvm_ind) # 선택된 변수중 실제 true 변수 비율

ssvm_f1score = 2 * (ssvm_precision * ssvm_recall) / (ssvm_precision + ssvm_recall) 


elssvm_ind = which(elssvm_result$theta != 0)
elssvm_select = intersect(true_ind, elssvm_ind) # true 변수중에 선택된거

elssvm_recall = length(elssvm_select)/length(true_ind) # 실제 true중 선택된 변수 비율
elssvm_precision = length(elssvm_select) / length(elssvm_ind)

elssvm_f1score = 2 * (elssvm_precision * elssvm_recall) / (elssvm_precision + elssvm_recall) 

###################################################################################

###########################################################################
##### 2022 - 03 - 00 ~ simulation, linear/nonlinear  data simulation ######
###########################################################################

# perform list
ssvm_er = c()
grid_er = c()
elssvm_er = c()

# select list
ssvm_f1score = c()
grid_f1score = c()
elssvm_f1score = c()

ssvm_recall = c()
grid_recall = c()
elssvm_recall = c()

# sub
opt_ind = c()

gamma = seq(0.1, 0.9, length.out = 5)




################################################################################


for (i in 1:10) {
cat("simulation_number:", i, "\n")

train_data = bsim_gen(150, 10, seed = i, type = "linearA")
train_x = scale(train_data$x)
train_y = as.numeric(as.character(train_data$y))

test_data = bsim_gen(10000, 10, seed = i, type = "linearA")
test_x = scale(test_data$x)
test_y = as.numeric(as.character(test_data$y))

ssvm_out = ssvm(x = train_x, y = train_y, nfolds = 5,lambda_seq = 2^seq(-10, 10, length.out = 100), 
                lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, cv_type = "original", isCombined = TRUE,
                kernel = "gaussian", kparam = 1,
                scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, verbose = 1, epsilon_D = 1e-4)

true_ind = which(train_data$true == 1)
ssvm_ind = which(ssvm_out$theta != 0)

ssvm_select = intersect(true_ind, ssvm_ind) # true 변수중에 선택된거

ssvm_recall[i] = length(ssvm_select) / length(true_ind) # 실제 true중 선택된 변수 비율
ssvm_precision = length(ssvm_select) / length(ssvm_ind) # 선택된 변수중 실제 true 변수 비율

ssvm_f1score[i] = 2 * (ssvm_precision * ssvm_recall[i]) / (ssvm_precision + ssvm_recall[i]) 

ssvm_pred = predict.ssvm(object = ssvm_out, newx = test_x)
ssvmresult = score(y_test = test_y, pred = ssvm_pred$class)

ssvm_er[i] = ssvmresult$error_rate

    for (j in 1:5){
            elssvm_out = elssvm(x = train_x, y = train_y, nfolds = 5,lambda_seq = 2^seq(-10, 10, length.out = 100), 
                                lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, cv_type = "original", isCombined = TRUE,
                                kernel = "gaussian", kparam = 1,
                                scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, verbose = 1, epsilon_D = 1e-4, gamma = gamma[j])
        
            elssvm_ind = which(elssvm_out$theta != 0)

            elssvm_select = intersect(true_ind, elssvm_ind) # true 변수중에 선택된거

            grid_elssvm_recall = length(elssvm_select) / length(true_ind) # 실제 true중 선택된 변수 비율
            elssvm_precision = length(elssvm_select) / length(elssvm_ind) # 선택된 변수중 실제 true 변수 비율
            grid_f1score[j] = 2 * (elssvm_precision * grid_elssvm_recall) / (elssvm_precision + grid_elssvm_recall) 
            grid_recall[j] = grid_elssvm_recall

            grid_elssvm_pred = predict.elssvm(object = elssvm_out, newx = test_x)
            grid_elssvmresult = score(y_test = test_y, pred = grid_elssvm_pred$class)
            grid_er[j] = grid_elssvmresult$error_rate
        }

opt_ind = which(grid_er == min(grid_er))[1]

elssvm_recall[i] = grid_recall[opt_ind]
elssvm_f1score[i] = grid_f1score[opt_ind]

elssvm_er[i] = grid_er[opt_ind]
}

################################################################################


# 성능비교
round(mean(ssvm_er), 4)
round(sd(ssvm_er) / sqrt(150), 4)
round(mean(elssvm_er), 4)
round(sd(elssvm_er) / sqrt(150), 4)

# 변수선택
round(mean(ssvm_f1score), 4)
round(sd(ssvm_f1score)/sqrt(150), 4)
round(mean(elssvm_f1score), 4)
round(sd(elssvm_f1score)/sqrt(150), 4)


################################################################################



