#%%
library(mlbench)
library(MASS)
#%%

# binary case, class 2 generate simulation

bsim_gen = function(n, p, seed = NULL, type = c("bayes", "linearA", "linearB", "neuralnet"))
{
    call = match.call()
    type = match.arg(type)

    if (!is.null(seed)) {
    set.seed(seed)
}
if (type == "bayes") {
    dat = mlbench.2dnormals(n = n, cl = 2, sd = 1)
    X_true = dat$x
    y = dat$classes
    r = ncol(X_true)
    X_noise = matrix(rnorm(n * (p - r), sd = 1), nrow = n, ncol = p - r)
    X = cbind(X_true, X_noise)

    out = list()
    out$x = X
    out$y = y
    out$true = rep(c(1, 0), c(r, p - r)) # 1이면 true값, 0이면 noise
}

if (type == "linearA") {
rho = matrix(0, 8, 8)
diag(rho) = 1

for (i in 1:8){
    for (j in 1:8){
        if (i == j){
            rho[i,j] = 1
        }
        else{
            rho[i,j] = 0.5^(abs(i - j))
        }
}
}

mutemp = rep(0, 8)
X = mvrnorm(n = n, mu = mutemp, Sigma = rho)
beta = c(3, 1.5, 1, 1, 2, 1, 1, 1)

for (i in 1:8){
    X[,i] = X[,i] * beta[i]
}

prob = exp(X) / (1 + exp(X))
y = rbinom(n = n, size = 1 ,prob = prob)

y_temp = factor(y)
classname = levels(y_temp)
n_class = length(classname)

y_int = integer(length(y))
for (j in 1:n_class) {y_int[which(y_temp %in% classname[j])] = j}
if (is(y, "numeric")) {classname = as.numeric(classname)}
y_int = ifelse(y_int == 1, 1, -1)

r = ncol(X)
X_noise = matrix(rnorm(n * (p - r), sd = 1), nrow = n, ncol = p - r)
X = cbind(X, X_noise)

out = list()

out$x = X
out$y = y
out$true = rep(c(1, 0), c(r, p - r)) # 1이면 true값, 0이면 noise

}

if (type == "linearB") {
rho = matrix(0, 10, 10)
diag(rho) = 1

for (i in 1:10){
    for (j in 1:10){
        if (i == j){
            rho[i,j] = 1
        }
        else{
            rho[i,j] = 0.5
        }
}
}

mutemp = rep(0, 10)
X = mvrnorm(n = n, mu = mutemp, Sigma = rho)
beta = c(2, 2, 2, 2, 2, 1, 1, 1, 1, 1)


for (i in 1:10){
    X[,i] = X[,i] * beta[i]
}

prob = exp(X) / (1 + exp(X))
y = rbinom(n = n, size = 1 ,prob = prob)

y_temp = factor(y)
classname = levels(y_temp)
n_class = length(classname)

y_int = integer(length(y))
for (j in 1:n_class) {y_int[which(y_temp %in% classname[j])] = j}
if (is(y, "numeric")) {classname = as.numeric(classname)}
y_int = ifelse(y_int == 1, 1, -1)

if (p > ncol(X))
{
r = ncol(X)
X_noise = matrix(rnorm(n * (p - r), sd = 1), nrow = n, ncol = p - r)
X = cbind(X, X_noise)
}
else {
    r = 10
}

out = list()

out$x = X
out$y = y
out$true = rep(c(1, 0), c(r, p - r)) # 1이면 true값, 0이면 noise

}


if (type == "neuralnet") {


r = 3
rho = matrix(0, r, r)


for (i in 1:3){
    for (j in 1:3){
        if (i == j){
            rho[i,j] = 1
        }
        else{
            rho[i,j] = 0.3^(abs(i - j))
        }
    }
    }


mutemp = rep(0, 3)

X_tmp = mvrnorm(n = n, mu = mutemp, Sigma = rho)
X_tmp = cbind(X_tmp, X_tmp[,1] * X_tmp[,2], X_tmp[,1] * X_tmp[,3], X_tmp[,2] * X_tmp[,3])

sigmoid = function(x) {
return(1 / (1 + exp(-x)))
}
    
    beta_mat1 = matrix(c(2.1, 1.4, -1.2, 2.1, 1.0, 2.4,
                        -1.1, 2.4, 1.5, 3.4, 2.5, 0.7),
                        nrow = 6, byrow = FALSE)
    
    node1_mat = drop(X_tmp %*% beta_mat1) + 1
    layer1 = matrix(sigmoid(node1_mat), nrow = n, ncol = 2)
    y = apply(layer1, 1, which.max) 

    noise = matrix(rnorm(n * (p - r), 0, 1.0), n, (p - r))
    X = cbind(X_tmp[,c(-4, -5, -6)], noise)
    true_vec = c(rep(1, r), rep(0, p - r))
    out = list()
    
    out$x = X
    out$y = y
    out$true = true_vec
    
}
return(out)
}

