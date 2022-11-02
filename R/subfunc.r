
# is used to compute a list of anova kernel matrices.
# main.kernel
kernelMatrix = function(x, y, kernel = "gaussian", kparam = 1.0) {
  
  x = as.matrix(x)
  y = as.matrix(y)
  p = ncol(x)
  
  if (NCOL(x) == 0) {
    x = matrix(0, nrow = nrow(x), ncol = 1)
  }
  
  if (NCOL(y) == 0) {
    y = matrix(0, nrow = nrow(y), ncol = 1)
  }
  
  if (kernel == "poly") {
    K = (x %*% t(y) + 1.0)^kparam
  } else if(kernel == "gaussian" | kernel == "gaussian2") {
    normx = rowSums(x^2)
    normy = rowSums(y^2)
    temp = x %*% t(y)
    temp = (-2.0 * temp) + outer(normx, rep(1.0, nrow(y)), "*") + outer(rep(1.0, nrow(x)), normy, "*")
    K = exp(-temp * kparam)
    # obj = kernelMatrix(rbfdot(sigma = kparam), x, y)
  } else if (kernel == "spline") {
    K = 0
    for (d in 1:p) {
      K_temp = spline_kernel(x[, d, drop = FALSE], y[, d, drop = FALSE])
      K = K + K_temp$K1 + K_temp$K2
    }
  } else if (kernel == "linear") {
    K = tcrossprod(x, y)
  } else if (kernel == "anova_gaussian") {
    K = 0
    for (d in 1:p) {
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = kernelMatrix(A, B, kernel = "gaussian", kparam = kparam)
      K = K + K_temp
    }
  } else {
    K = NULL
  }
  return(K)
}



spline_kernel = function(x, y) 
{
  x = as.matrix(x)
  y = as.matrix(y)
  K1x = (x - 1 / 2)
  K1y = (y - 1 / 2)
  K2x = (K1x^2 - 1 / 12) / 2
  K2y = (K1y^2 - 1 / 12) / 2
  ax = x %x% matrix(1, 1, nrow(y)) 
  ay = y %x% matrix(1, 1, nrow(x))
  b = abs(ax - t(ay))
  K1 = K1x %x% t(K1y)
  K2 = K2x %x% t(K2y) - ((b - 1 / 2)^4 - (b - 1 / 2)^2 / 2 + 7 / 240) / 24
  list(K1 = K1, K2 = K2)
}


# kernel$type: t 붙은건 combined method, 그냥은 separated method
# 'spline' for additive models with linear and smooth parts separated.
# 'spline-t' for additive models with linear and smooth parts combined.
# 'spline2' for two-way interaction models with linear and smooth parts separated.
# 'spline-t2' for two-way interaction models with linear and smooth parts combined.
make_anovaKernel = function(x, y, kernel, kparam)
{
  x = as.matrix(x)
  y = as.matrix(y)
  dimx = ncol(x)
  
  # calculate anova kernels for main effects
  if (kernel == "spline") {
    # assign the number of anova kernels
    numK = 2 * dimx
    # list of kernel matrices
    anova_kernel = vector(mode = "list", numK)
    # list of kernel coordinate indices
    kernelCoord = vector(mode = "list", numK)
    index = 0
    
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep="")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep="")
    } 
    
  } else if (kernel == 'spline2') {
    numK = (2 * dimx) + (2 * dimx * (2 * dimx - 1) / 2 - dimx)
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    # main effects
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep = "")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep = "")
    }  
    # two-way interactions
    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A_linear = as.matrix(anova_kernel[[2 * i - 1]])
        A_smooth = as.matrix(anova_kernel[[2 * i]])
        B_linear = as.matrix(anova_kernel[[2 * j - 1]])
        B_smooth = as.matrix(anova_kernel[[2 * j]])            
        anova_kernel[[index]] = A_linear * B_linear
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " linear", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_linear * B_smooth
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " smooth", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_smooth * B_linear
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " linear", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_smooth * B_smooth
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " smooth", sep = "")
      }
    }
  } else if (kernel == "spline-t") {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      kernelCoord[[index]] = paste("x", d, sep = "")
    }
  } else if (kernel == 'spline-t2') {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      kernelCoord[[index]] = paste("x", d, sep = "")
    }
    
    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A = anova_kernel[[i]]
        B = anova_kernel[[j]]
        anova_kernel[[index]] = A * B
        kernelCoord[[index]] = paste("x", i, " x", j, sep = "")
      }
    }
  } else if (kernel == "gaussian2") {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      anova_kernel[[index]] = kernelMatrix(A, B, kernel, kparam)
      kernelCoord[[index]] = paste("x", d, sep = "")
    }
    
    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A = anova_kernel[[i]]
        B = anova_kernel[[j]]
        anova_kernel[[index]] = A * B
        kernelCoord[[index]] = paste("x", i, " x", j, sep = "")
      }
    }
  } else {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    for (d in 1:dimx) {
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      anova_kernel[[d]] = kernelMatrix(A, B, kernel, kparam)
      kernelCoord[[d]] = paste("x", d, sep = "")
    }
  }
  return(list(K = anova_kernel, coord = kernelCoord, numK = numK, kernel = kernel, kparam = kparam))
}

# is used to combine anova kernel matrices with weights determined by theta values. The default theta vector is the vector of ones.
combine_kernel = function(anova_kernel, theta = rep(1, anova_kernel$numK))
{
  K = 0
  for (d in 1:anova_kernel$numK) {
    K = (K + theta[d] * anova_kernel$K[[d]])
  }
  return(K)
}


data_split = function(y, nfolds, seed = length(y))
{
  # k: the number of classes
  y = as.factor(y)
  n_data = length(y)
  n_class = length(levels(y))
  class_size = table(y)
  classname = names(class_size)
  ran = rep(0, n_data) 
  if ((min(class_size) < nfolds) & (nfolds != n_data))
  {
    warning('The given fold is bigger than the smallest class size. \n Only a fold size smaller than the minimum class size \n or the same as the sample size (LOOCV) is supported.\n')
    return(NULL)
  }
  
  if (min(class_size) >= nfolds) {
    set.seed(seed)
    for (j in 1:n_class) {  
      ran[y == classname[j]] = ceiling(sample(class_size[j]) / (class_size[j] + 1) * nfolds) 
    }
  }
  else if (nfolds == n_data) {
    ran = 1:n_data
  }
  return(ran)
}



strong_heredity = function(main_effect, interaction)
{
  if (sum(interaction) != 0) {
    p = length(main_effect)
    comb_mat = combn(1:p, 2)
    ind_mat = comb_mat[, interaction == 1, drop = FALSE]
    for (i in 1:nrow(ind_mat)) {
      ind = ind_mat[i, ]
      main_effect[ind] = 1
    }
  }
  res = c(main_effect, interaction)
  return(res)
}


# interaction_graph = function(comb, p, min = 3)
# {
#   int_mat = Matrix::Matrix(0, nrow = p, ncol = p)
#   int_mat[t(comb)] = 1
#   int_mat = Matrix::t(int_mat) + int_mat
#   g = graph_from_adjacency_matrix(int_mat, mode = "undirected")
#   cliques_list = max_cliques(g, min = 3)
#   return(cliques_list)
# }

interaction_kernel = function(x, u, kernel, kparam, active_set, interaction_set)
{
  if (!is.matrix(x)) {
    x = as.matrix(x)
  } else {
    x = as.matrix(x)
  }
  u = as.matrix(u)
  dimx = ncol(x)
  
  scaled_kernel = function(x, u, kernel, kparam, active_set, index)
  {
    X1 = matrix(rowSums(x[, active_set, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u))
    U1 = matrix(rowSums(u[, active_set, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u), byrow = TRUE)
    X2 = matrix(rowSums(x[, index, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u))
    U2 = matrix(rowSums(u[, index, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u), byrow = TRUE)
    K = exp(-kparam * ((X1 + U1) - (X2 + U2)))
    K1 = exp(-kparam * (X1 + U1))
    K2 = exp(-kparam * (X2 + U2))
    K_mat = kernelMatrix(x[, index, drop = FALSE], u[, index, drop = FALSE], kernel = kernel, kparam = kparam)
    res = K * K_mat - K1
    return(list(res = res, K = K, K1 = K1, K2 = K2))
  }
  
  main_effects = vector(mode = "list", dimx)
  
  for (j in 1:length(active_set)) {
    temp_kernel = scaled_kernel(x, u, kernel = kernel, kparam = kparam, active_set = active_set, index = active_set[j])
    main_effects[[active_set[j]]] = temp_kernel$res
  }
  
  interaction_kernel = 0
  if (ncol(interaction_set) != 0) {
    
    for (i in 1:ncol(interaction_set)) {
      ind = interaction_set[, i]
      interaction_kernel = interaction_kernel + ((main_effects[[ind[1]]]) * (main_effects[[ind[2]]])) / temp_kernel$K1
    }
  }
  
  K = temp_kernel$K1 + Reduce("+", main_effects[active_set]) + interaction_kernel
  
  return(K)
}

# is used to find nonzero entries of a given matrix and their positions in the matrix for a compact representation of elements of QP problem.
find_nonzero = function(Amat)
{
  nr = nrow(Amat)
  nc = ncol(Amat)
  Amat_compact = matrix(0, nr, nc)
  Aind = matrix(0, nr + 1, nc)
  for (j in 1:nc) {
    index = (1:nr)[Amat[, j] != 0]
    number = length(index)
    Amat_compact[1:number, j] = Amat[index, j]
    Aind[1, j] = number
    Aind[2:(number+1), j] = index
  }
  max_number = max(Aind[1, ])
  Amat_compact = Amat_compact[1:max_number, ]
  Aind = Aind[1:(max_number + 1), ]
  return(list(Amat_compact = Amat_compact, Aind = Aind))
}

fixit = function(A, epsilon = .Machine$double.eps, is_diag = FALSE)
{
  if (is_diag) {
    d = diag(A)
    tol = epsilon
    eps = max(tol * max(d), 0)
    d[d < eps] = eps
    Q = diag(d)
  } else {
    eig = eigen(A, symmetric = TRUE)
    tol = epsilon
    eps = max(tol * abs(eig$values[1]), 0)
    eig$values[eig$values < eps] = eps
    Q = eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  }
  return(Q)
}



unit.scale <- function(x, limit = NULL) # spline kernel 사용시 x, data scale
{
  x <- as.matrix(x)   
  if (is.null(limit))
  {
    r <- apply(x,2,range)
    x <- scale(x, r[1,] - .05*(r[2,]-r[1,]), 1.1*(r[2,]-r[1,]))
  }
  else
  {
    r <- limit 
    x <- scale(x, r[1,] - .05*(r[2,]-r[1,]), 1.1*(r[2,]-r[1,]))
    x[x<0] <- 0
    x[x>1] <- 1
  }
  return(x)
}


score = function(pred, y_test){

confusion = confusionMatrix(as.factor(pred), as.factor(test_y))

accuracy = confusion$overall["Accuracy"] # 정분류율
error_rate = 1 - accuracy

cat("\naccuracy:",accuracy, "\nerror_rate:",error_rate)
  return(list(accuracy = accuracy, error_rate = error_rate))
}

###
# install.packages("ggplot")
# library(ggplot2)


# prost_mean = c(0.0393, 0.0303, 0.0394, 0.0255)
# parkin_mean = c(0.1991, 0.1675, 0.2569, 0.1964)
# onar_mean = c(0.2510, 0.2294, 0.1684, 0.1471)
# iono_mean = c(0.1333, 0.1198, 0.0819, 0.0768)
# forest_mean = c(0.0268, 0.0186, 0.0277, 0.0189)
# heart_mean = c(0.1893, 0.1713, 0.1676, 0.1498)

# mean_df = data.frame(cbind(prost_mean, parkin_mean, sonar_mean, iono_mean, forest_mean, heart_mean))

# olnames(mean_df) = c("Prostate Cancer" , "Parkinson", "Sonar", "Ionosphere", "Forest Fires", "Heart Clinical")

# ggplot(data=mean_df, aes(y="miclassification rate")) +
#    geom_bar(stat="identity", position=position_dodge())
















