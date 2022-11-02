# EL-COSSO-SVM
Variable selection for nonlinear support vector machines via elastic net penalty

## 1. INSTALLATIO

## 2. EXAMPLE

```{r}
# Generation of simulated data with linear decision boundary
> require(elssvm)
> n = 100; p = 10; 
> data = elssvm:::sim_gen(n = n, p = p, type = "linear")
> x = scale(data$x)
> y = data$y
> sigma = kernlab::sigest(y ~ x, scaled = FALSE)[3]

```