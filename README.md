# PCS

## Partial correlation screening approach for controlling the false positive rate in sparse Gaussian Graphical Models

The R package contains functions used in the following [paper](https://www.nature.com/articles/s41598-019-53795-x):

Lafit G, Tuerlinckx F, Myin-Germeys I, Ceulemans E. A Partial Correlation Screening Approach for Controlling the False Positive Rate in Sparse Gaussian Graphical Models. Sci Rep. 2019;9(1):17759. Published 2019 Nov 28. doi:10.1038/s41598-019-53795-x

## Installation
You can install the development version from GitHub with:

```
library(devtools)
devtools::install_github("ginettelafit/PCS", force = T)
```
## Example - simulate data 

To illustrate how the PCS approach works, we simulated data by drawing 100 independent observations from a multivariate Gaussian distribution with mean zero and partial correlation matrix ![](http://latex.codecogs.com/svg.latex?%5Cboldsymbol%7B%5CGamma%7D). We consider the dimension to be euqual to 20. 
We assumed 2 neighbor Chain Graph for the population partial correlation matrix, in which ![](http://latex.codecogs.com/svg.latex?%5Crho_%7Bii%7CV%5Csetminus%20%5C%7Bi%5C%7D%7D%3D1) and ![](http://latex.codecogs.com/svg.latex?%5Csetminus%20%5C%7Bi%2Ci&plus;1%5C%7D%7D%3D%5Crho_%7Bi-1%2Ci%7CV%5Csetminus%20%5C%7Bi%2Ci-1%5C%7D%7D%3D-.4), and all other edges are set to 0.

```
library(Rlab)
library(psych)
library(MASS)
library(huge)

p = 20 # Set dimenensionality
n = 100 # Set sample size
 
# Generate the covariance (Sigma) matrix & simulate data N(0,Sigma) 

# Model 1: : 2 neighbor Chain Graph

Omega = diag(p)

for (i in 1:p){
Omega[i,i-1] = 0.4
Omega[i-1,i] = 0.4
}

Sigma = solve(Omega)

R = -cov2cor(Omega)
diag(R) = 1

save(R, file = "pcor_1_p_20.RData")

# Simulate Data
x = mvrnorm(n, rep(0, p), Sigma)
x = scale(x, center = TRUE, scale = F)
save(x, file = "sample_1_p_20_n_100.RData")
```

## Estimate GGM using Glasso, SPACE, and nodewise regression using lasso

```
# Glasso 10-folds CV1

  load(file = "sample_1_p_20_n_100.RData")
  rho.list = cv_glasso(x,fold=10)    
  rho = unlist(rho.list$rho_opt)
  list.glasso.cv = as.matrix(huge(x,lambda = rho,scr = F,method = "glasso")$icov[[1]])
  lambda.opt = rho
  save(lambda.opt, file = "lambda_glasso_cv_1_p_20_n_100.RData")
  save(list.glasso.cv, file = "glasso_cv_1_p_20_n_100.RData")

# Glasso 10-folds CV1-1se

  rho.ls1 = unlist(rho.list$rho_ls1)
  list.glasso.cv.ls1 = as.matrix(huge(x,lambda = rho.ls1,scr = F,method = "glasso")$icov[[1]])
  lambda.opt = rho.ls1
  save(lambda.opt, file = "lambda_glasso_cv_ls1_1_p_20_n_100.RData")
  save(list.glasso.cv.ls1, file = "glasso_cv_ls1_1_p_20_n_100.RData")

# Glasso 10-folds CV2

  rho = unlist(rho.list$rho2_opt)
  list.glasso.cv2 = as.matrix(huge(x,lambda = rho,scr = F,method = "glasso")$icov[[1]])
  lambda.opt = rho
  save(lambda.opt, file = "lambda_glasso_cv2_1_p_20_n_100.RData")
  save(list.glasso.cv2, file = "glasso_cv2_1_p_20_n_100.RData")

# Glasso 10-folds CV2-1se

  rho.ls1 = unlist(rho.list$rho2_ls1)
  list.glasso.cv2.ls1 = as.matrix(huge(x,lambda = rho.ls1,scr = F,method = "glasso")$icov[[1]])
  lambda.opt = rho.ls1
  save(lambda.opt, file = "lambda_glasso_cv2_ls1_1_p_20_n_100.RData")
  save(list.glasso.cv2.ls1, file = "glasso_cv2_ls1_1_p_20_n_100.RData")

# Glasso EBIC
# Glasso + EBIC estimated using the function huge.select from the R package huge
library(huge)

  load(file = "sample_1_p_20_n_100.RData")
  fit.ebic = huge.select(huge(x,lambda=seq(0.001,max(abs(cov(x))),length=100),method = "glasso"),
  ebic.gamma = 0.5,criterion = "ebic")
  lambda.opt = fit.ebic$opt.lambda
  list.glasso.ebic = as.matrix(fit.ebic$opt.icov)
  save(lambda.opt, file = "lambda_glasso_ebic_1_p_20_n_100.RData")
  save(list.glasso.ebic, file = "glasso_ebic_1_p_20_n_100.RData")

# Glasso BIC
# Glasso + BIC estimated using the function huge.select from the R package huge
library(huge)

  load(file = "sample_1_p_20_n_100.RData")
  fit.bic = huge.select(huge(x,lambda=seq(0.001,max(abs(cov(x))),length=100),method = "glasso"),
  ebic.gamma = 0,criterion = "ebic")
  lambda.opt = fit.bic$opt.lambda
  list.glasso.bic = as.matrix(fit.bic$opt.icov)
  save(lambda.opt, file = "lambda_glasso_bic_1_p_20_n_100.RData")
  save(list.glasso.bic, file = "glasso_bic_1_p_20_n_100.RData")

# SPACE with 10-folds CV

  load(file = "sample_1_p_20_n_100.RData")
  lambda.list = cv_space(x,fold=10)
  lambda.opt = lambda.list$alpha_opt
  R.opt = Psi_space(x,lambda.opt)
  list.space.cv = PCS_GGM(x,R.opt,0)
  save(lambda.opt, file = "lambda_space_cv_1_p_20_n_100.RData")
  save(list.space.cv, file = "space_cv_1_p_20_n_100.RData")

  lambda.ls1 = lambda.list$alpha_ls1
  R.opt = Psi_space(x,lambda.ls1)
  list.space.cv.ls1 = PCS_GGM(x,R.opt,0)
  lambda.opt = lambda.ls1
  save(lambda.opt, file = "lambda_space_cv_ls1_1_p_20_n_100.RData")
  save(list.space.cv.ls1, file = "space_cv_ls1_1_p_20_n_100.RData")

# SPACE with finite sample results result

  load(file = "sample_1_p_20_n_100.RData")
  lambda.opt = sqrt(nrow(x))*qnorm(1-0.05/(2*ncol(x)^2))
  R.opt = Psi_space(x,lambda.opt)
  list.space.alpha = PCS_GGM(x,R.opt,0)
  save(lambda.opt, file = "lambda_space_alpha_1_p_20_n_100.RData")
  save(list.space.alpha, file = "space_alpha_1_p_20_n_100.RData")

# SPACE with BIC

  load(file = "sample_1_p_20_n_100.RData")
  lambda.opt = BIC_space(x)
  R.opt = Psi_space(x,lambda.opt)
  list.space.bic = PCS_GGM(x,R.opt,0)
  save(lambda.opt, file = "lambda_space_bic_1_p_20_n_100.RData")
  save(list.space.bic, file = "space_bic_1_p_20_n_100.RData")

# Nodewise regression 'AND' with 10-folds CV

  load(file = "sample_1_p_20_n_100.RData")
  lambda.fit = cv_Neigh(x,10)
  lambda.opt = lambda.fit[[1]]
  list.nei.cv = neigh_and_lambda(x,lambda.opt)
  save(lambda.opt, file = "lambda_nei_and_cv_1_p_20_n_100.RData")
  save(list.nei.cv, file = "nei_and_cv_1_p_20_n_100.RData")

# Nodewise regression 'AND' with 10-folds CV-ls1

  list.nei.cv = neigh_or_lambda(x,lambda.opt)
  save(lambda.opt, file = "lambda_nei_or_cv_1_p_20_n_100.RData")
  save(list.nei.cv, file = "nei_or_cv_1_p_20_n_100.RData")

# Nodewise regression 'OR' with 10-folds CV

  lambda.ls1 = lambda.fit[[2]]
  list.nei.cv.ls1 = neigh_and_lambda(x,lambda.ls1)
  lambda.opt = lambda.ls1
  save(lambda.opt, file = "lambda_nei_and_cv_ls1_1_p_20_n_100.RData")
  save(list.nei.cv.ls1, file = "nei_and_cv_ls1_1_p_20_n_100.RData")

# Nodewise regression 'OR' with 10-folds CV-ls1

  list.nei.cv.ls1 = neigh_or_lambda(x,lambda.ls1)
  save(lambda.opt, file = "lambda_nei_or_cv_ls1_1_p_20_n_100.RData")
  save(list.nei.cv.ls1, file = "nei_or_cv_ls1_1_p_20_n_100.RData")

# Nodewise regression 'AND' with finite sample result

  load(file = "sample_1_p_20_n_100.RData")
  lambda.opt = rep((1/(sqrt(nrow(x))))*qnorm(1-0.05/(2*ncol(x)^2)),ncol(x))
  list.nei.alpha = neigh_and_lambda(x,lambda.opt)
  save(lambda.opt, file = "lambda_nei_and_alpha_1_p_20_n_100.RData")
  save(list.nei.alpha, file = "nei_and_alpha_1_p_20_n_100.RData")

# Nodewise regression 'OR' with finite sample result

  list.nei.alpha = neigh_or_lambda(x,lambda.opt)
  save(lambda.opt, file = "lambda_nei_or_alpha_1_p_20_n_100.RData")
  save(list.nei.alpha, file = "nei_or_alpha_1_p_20_n_100.RData")

# Nodewise regression 'AND' with BIC

  load(file = "sample_1_p_20_n_100.RData")
  lambda.opt = BIC_Neigh(x)
  list.nei.bic = neigh_and_lambda(x,lambda.opt)
  save(lambda.opt, file = "lambda_nei_and_bic_1_p_20_n_100.RData")
  save(list.nei.bic, file = "nei_and_bic_1_p_20_n_100.RData")

# Nodewise regression 'OR' with BIC

  list.nei.bic = neigh_or_lambda(x,lambda.opt)
  save(lambda.opt, file = "lambda_nei_or_bic_1_p_20_n_100.RData")
  save(list.nei.bic, file = "nei_or_bic_1_p_20_n_100.RData")

# Ridge regression with 10-folds CV

  load(file = "sample_1_p_20_n_100.RData")
  lambda.opt = cv_ridge(x,10)
  list.ridge.cv = ridge_lambda(x,lambda.opt)
  save(lambda.opt, file = "lambda_ridge_cv_1_p_20_n_100.RData")
  save(list.ridge.cv, file = "ridge_cv_1_p_20_n_100.RData")
```

## Control the false positive rate using PCS

```
# Glasso 10-folds CV1

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "glasso_cv_1_p_20_n_100.RData")
  load(file = "lambda_glasso_cv_1_p_20_n_100.RData")
  R.opt = -cov2cor(list.glasso.cv)
  diag(R.opt) = 1
  tau.opt = cv_PCS_glasso(x,lambda.opt,fold=10)
  list.glasso.cv.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_glasso_cv_1_p_20_n_100.RData")
  save(list.glasso.cv.pcs, file = "PCS_glasso_cv_1_p_20_n_100.RData")

# Glasso 10-folds CV1-1se

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "glasso_cv_ls1_1_p_20_n_100.RData")
  load(file = "lambda_glasso_cv_ls1_1_p_20_n_100.RData")
  R.opt = -cov2cor(list.glasso.cv.ls1)
  diag(R.opt) = 1
  tau.opt = cv_PCS_glasso(x,lambda.opt,fold=10)
  list.glasso.cv.ls1.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_glasso_cv_ls1_1_p_20_n_100.RData")
  save(list.glasso.cv.ls1.pcs, file = "PCS_glasso_cv_ls1_1_p_20_n_100.RData")

# Glasso 10-folds CV2

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "glasso_cv2_1_p_20_n_100.RData")
  load(file = "lambda_glasso_cv2_1_p_20_n_100.RData")
  R.opt = -cov2cor(list.glasso.cv2)
  diag(R.opt) = 1
  tau.opt = cv_PCS_glasso(x,lambda.opt,fold=10)
  list.glasso.cv2.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_glasso_cv2_1_p_20_n_100.RData")
  save(list.glasso.cv2.pcs, file = "PCS_glasso_cv2_1_p_20_n_100.RData")

# Glasso 10-folds CV2-1se

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "glasso_cv2_ls1_1_p_20_n_100.RData")
  load(file = "lambda_glasso_cv2_ls1_1_p_20_n_100.RData")
  R.opt = -cov2cor(list.glasso.cv2.ls1)
  diag(R.opt) = 1
  tau.opt = cv_PCS_glasso(x,lambda.opt,fold=10)
  list.glasso.cv2.ls1.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_glasso_cv2_ls1_1_p_20_n_100.RData")
  save(list.glasso.cv2.ls1.pcs, file = "PCS_glasso_cv2_ls1_1_p_20_n_100.RData")

# Glasso EBIC

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "glasso_ebic_1_p_20_n_100.RData")
  load(file = "lambda_glasso_ebic_1_p_20_n_100.RData")
  R.opt = -cov2cor(list.glasso.ebic)
  diag(R.opt) = 1
  tau.opt = cv_PCS_glasso(x,lambda.opt,fold=10)
  list.glasso.ebic.pcs = PCS_GGM(x,R.opt,tau.opt)
  save(tau.opt, file = "tau_glasso_ebic_1_p_20_n_100.RData")
  save(list.glasso.ebic.pcs, file = "PCS_glasso_ebic_1_p_20_n_100.RData")

# Glasso BIC

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "glasso_bic_1_p_20_n_100.RData")
  load(file = "lambda_glasso_bic_1_p_20_n_100.RData")
  R.opt = -cov2cor(list.glasso.bic)
  diag(R.opt) = 1
  tau.opt = cv_PCS_glasso(x,lambda.opt,fold=10)
  list.glasso.bic.pcs = PCS_GGM(x,R.opt,tau.opt)
  save(tau.opt, file = "tau_glasso_bic_1_p_20_n_100.RData")
  save(list.glasso.bic.pcs, file = "PCS_glasso_bic_1_p_20_n_100.RData")

# SPACE with 10-folds CV

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "space_cv_1_p_20_n_100.RData")
  load(file = "lambda_space_cv_1_p_20_n_100.RData")
  R.opt = Psi_space(x,lambda.opt)
  #R.opt = -cov2cor(list.space.cv[[2]])
  #diag(R.opt) = 1
  tau.opt = cv_PCS_space(x,lambda.opt,10)
  list.space.cv.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_space_cv_1_p_20_n_100.RData")
  save(list.space.cv.pcs, file = "PCS_space_cv_1_p_20_n_100.RData")

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "space_cv_ls1_1_p_20_n_100.RData")
  load(file = "lambda_space_cv_ls1_1_p_20_n_100.RData")
  R.opt = Psi_space(x,lambda.opt)
  #R.opt = -cov2cor(list.space.cv.ls1[[2]])
  #diag(R.opt) = 1
  tau.opt = cv_PCS_space(x,lambda.opt,10)
  list.space.cv.ls1.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_space_cv_ls1_1_p_20_n_100.RData")
  save(list.space.cv.ls1.pcs, file = "PCS_space_cv_ls1_1_p_20_n_100.RData")

# SPACE with finite sample results result

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "space_alpha_1_p_20_n_100.RData")
  load(file = "lambda_space_alpha_1_p_20_n_100.RData")
  R.opt = Psi_space(x,lambda.opt)
  #R.opt = -cov2cor(list.space.alpha[[2]])
  #diag(R.opt) = 1
  tau.opt = cv_PCS_space(x,lambda.opt,10)
  list.space.alpha.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_space_alpha_1_p_20_n_100.RData")
  save(list.space.alpha.pcs, file = "PCS_space_alpha_1_p_20_n_100.RData")

# SPACE with BIC

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "space_bic_1_p_20_n_100.RData")
  load(file = "lambda_space_bic_1_p_20_n_100.RData")
  R.opt = Psi_space(x,lambda.opt)
  #R.opt = -cov2cor(list.space.bic[[2]])
  #diag(R.opt) = 1
  tau.opt = cv_PCS_space(x,lambda.opt,10)
  list.space.bic.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_space_bic_1_p_20_n_100.RData")
  save(list.space.bic.pcs, file = "PCS_space_bic_1_p_20_n_100.RData")

# Nodewise regression 'AND' with 10-folds CV

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "nei_and_cv_1_p_20_n_100.RData")
  load(file = "lambda_nei_and_cv_1_p_20_n_100.RData")
  R.opt = -cov2cor(list.nei.cv[[2]])
  diag(R.opt) = 1
  tau.opt = cv_PCS_nei_and(x,lambda.opt,10) 
  list.nei.cv.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_nei_and_cv_1_p_20_n_100.RData")
  save(list.nei.cv.pcs, file = "PCS_nei_and_cv_1_p_20_n_100.RData")

# Nodewise regression 'AND' with 10-folds CV-ls1

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "nei_or_cv_1_p_20_n_100.RData")
  load(file = "lambda_nei_or_cv_1_p_20_n_100.RData")
  R.opt = -cov2cor(list.nei.cv[[2]])
  diag(R.opt) = 1
  tau.opt = cv_PCS_nei_or(x,lambda.opt,10)
  list.nei.cv.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_nei_or_cv_1_p_20_n_100.RData")
  save(list.nei.cv.pcs, file = "PCS_nei_or_cv_1_p_20_n_100.RData")

# Nodewise regression 'OR' with 10-folds CV

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "nei_and_cv_ls1_1_p_20_n_100.RData")
  load(file = "lambda_nei_and_cv_ls1_1_p_20_n_100.RData")
  R.opt = -cov2cor(list.nei.cv.ls1[[2]])
  diag(R.opt) = 1
  tau.opt = cv_PCS_nei_and(x,lambda.opt,10) 
  list.nei.cv.ls1.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_nei_and_cv_ls1_1_p_20_n_100.RData")
  save(list.nei.cv.ls1.pcs, file = "PCS_nei_and_cv_ls1_1_p_20_n_100.RData")

# Nodewise regression 'OR' with 10-folds CV-ls1

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "nei_or_cv_ls1_1_p_20_n_100.RData")
  load(file = "lambda_nei_or_cv_ls1_1_p_20_n_100.RData")
  R.opt = -cov2cor(list.nei.cv.ls1[[2]])
  diag(R.opt) = 1
  tau.opt = cv_PCS_nei_or(x,lambda.opt,10) 
  list.nei.cv.ls1.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_nei_or_cv_ls1_1_p_20_n_100.RData")
  save(list.nei.cv.ls1.pcs, file = "PCS_nei_or_cv_ls1_1_p_20_n_100.RData")

# Nodewise regression 'AND' with finite sample result

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "nei_and_alpha_1_p_20_n_100.RData")
  load(file = "lambda_nei_and_alpha_1_p_20_n_100.RData")
  R.opt = -cov2cor(list.nei.alpha[[2]])
  diag(R.opt) = 1
  tau.opt = cv_PCS_nei_and(x,lambda.opt,10)
  list.nei.alpha.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_nei_and_alpha_1_p_20_n_100.RData")
  save(list.nei.alpha.pcs, file = "PCS_nei_and_alpha_1_p_20_n_100.RData")

# Nodewise regression 'OR' with finite sample result

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "nei_or_alpha_1_p_20_n_100.RData")
  load(file = "lambda_nei_or_alpha_1_p_20_n_100.RData")
  R.opt = -cov2cor(list.nei.alpha[[2]])
  diag(R.opt) = 1
  tau.opt = cv_PCS_nei_or(x,lambda.opt,10)
  list.nei.alpha.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_nei_or_alpha_1_p_20_n_100.RData")
  save(list.nei.alpha.pcs, file = "PCS_nei_or_alpha_1_p_20_n_100.RData")

# Nodewise regression 'AND' with BIC

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "nei_and_bic_1_p_20_n_100.RData")
  load(file = "lambda_nei_and_bic_1_p_20_n_100.RData")
  R.opt = -cov2cor(list.nei.bic[[2]])
  diag(R.opt) = 1
  tau.opt = cv_PCS_nei_and(x,lambda.opt,10)
  list.nei.bic.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_nei_and_bic_1_p_20_n_100.RData")
  save(list.nei.bic.pcs, file = "PCS_nei_and_bic_1_p_20_n_100.RData")

# Nodewise regression 'OR' with BIC

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "nei_or_bic_1_p_20_n_100.RData")
  load(file = "lambda_nei_or_bic_1_p_20_n_100.RData")
  R.opt = -cov2cor(list.nei.bic[[2]])
  diag(R.opt) = 1
  tau.opt = cv_PCS_nei_or(x,lambda.opt,10)
  list.nei.bic.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_nei_or_bic_1_p_20_n_100.RData")
  save(list.nei.bic.pcs, file = "PCS_nei_or_bic_1_p_20_n_100.RData")

# Ridge regression with 10-folds CV

  load(file = "sample_1_p_20_n_100.RData")
  load(file = "ridge_cv_1_p_20_n_100.RData")
  load(file = "lambda_ridge_cv_1_p_20_n_100.RData")
  R.opt = -cov2cor(list.ridge.cv[[2]])
  diag(R.opt) = 1
  tau.opt = cv_PCS_ridge(x,lambda.opt,10)
  list.ridge.cv.pcs = PCS_GGM(x,R.opt,tau.opt) 
  save(tau.opt, file = "tau_ridge_cv_1_p_20_n_100.RData")
  save(list.ridge.cv.pcs, file = "PCS_ridge_cv_1_p_20_n_100.RData")
```
## Performance measures

The following plots show the heatmaps of the frequency with which the edges for the simulated data are set to zero by the different methods. White indicates that an edge was excluded from the network in all replications, whereas black reflects that the edge was always retained in the network.


The next heatmaps show the frequency with which the edges for the simulated data are set to zero by the different methods after applying PCS. 


The following plots show the networks estimated by the different methods. 


The next plots show the networks after applying PCS. 


