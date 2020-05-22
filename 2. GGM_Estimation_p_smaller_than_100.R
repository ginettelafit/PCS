############################################################
############################################################
#####         R Code by Ginette Lafit        ###############
#####         ginette.lafit@kuleuven.be      ###############
############################################################
############################################################

ls()

rm(list=ls())

set.seed(123)

############################################################
##### Load functions
############################################################

source(file="BIC_Neigh.R")
source(file="BIC_space.R")
source(file="cv_glasso.R")
source(file="cv_Neigh.R")
source(file="cv_space.R")
source(file="cv_ridge.R")
source(file="cv_PCS_glasso.R")
source(file="cv_PCS_nei_and.R")
source(file="cv_PCS_nei_or.R")
source(file="cv_PCS_ridge.R")
source(file="cv_PCS_space.R")
source(file="neigh_and_lambda.R")
source(file="neigh_or_lambda.R")
source(file="Psi_Screen_Beta.R")
source(file="Psi_space.R")
source(file="ridge_lambda.R")
source(file="PCS_GGM.R")

###########################################################
###########################################################
###########################################################
##### GGM Estimation
###########################################################
###########################################################
###########################################################

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

######################################################
######################################################
######################################################
######################################################

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

######################################################
######################################################
######################################################
######################################################

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

######################################################
######################################################
######################################################
######################################################

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

######################################################
######################################################
######################################################
######################################################

# SPACE with finite sample results result

  load(file = "sample_1_p_20_n_100.RData")
  lambda.opt = sqrt(nrow(x))*qnorm(1-0.05/(2*ncol(x)^2))
  R.opt = Psi_space(x,lambda.opt)
  list.space.alpha = PCS_GGM(x,R.opt,0)
  save(lambda.opt, file = "lambda_space_alpha_1_p_20_n_100.RData")
  save(list.space.alpha, file = "space_alpha_1_p_20_n_100.RData")

######################################################
######################################################
######################################################
######################################################

# SPACE with BIC

  load(file = "sample_1_p_20_n_100.RData")
  lambda.opt = BIC_space(x)
  R.opt = Psi_space(x,lambda.opt)
  list.space.bic = PCS_GGM(x,R.opt,0)
  save(lambda.opt, file = "lambda_space_bic_1_p_20_n_100.RData")
  save(list.space.bic, file = "space_bic_1_p_20_n_100.RData")

######################################################
######################################################
######################################################
######################################################

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

######################################################
######################################################
######################################################
######################################################

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

######################################################
######################################################
######################################################
######################################################

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

######################################################
######################################################
######################################################
######################################################

# Ridge regression with 10-folds CV

  load(file = "sample_1_p_20_n_100.RData")
  lambda.opt = cv_ridge(x,10)
  list.ridge.cv = ridge_lambda(x,lambda.opt)
  save(lambda.opt, file = "lambda_ridge_cv_1_p_20_n_100.RData")
  save(list.ridge.cv, file = "ridge_cv_1_p_20_n_100.RData")

######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################

###########################################################
###########################################################
###########################################################
##### GGM Estimation with PCS
###########################################################
###########################################################
###########################################################

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


######################################################
######################################################
######################################################
######################################################

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


######################################################
######################################################
######################################################
######################################################

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

######################################################
######################################################
######################################################
######################################################

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

######################################################
######################################################
######################################################
######################################################

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

######################################################
######################################################
######################################################
######################################################

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

######################################################
######################################################
######################################################
######################################################

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

######################################################
######################################################
######################################################
######################################################

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

######################################################
######################################################
######################################################
######################################################

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

######################################################
######################################################
######################################################
######################################################

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

######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################

