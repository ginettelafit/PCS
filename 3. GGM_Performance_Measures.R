############################################################
############################################################
#####         R Code by Ginette Lafit        ###############
#####         ginette.lafit@kuleuven.be      ###############
############################################################
############################################################

ls()

rm(list=ls())

library(psych)

library(corrplot)

library(qgraph)

library(gridExtra)

############################################################
############################################################
######          Classification Performance           #######
############################################################
############################################################

## Performance is a function that computes
# sensitivity and specificity of one matrix

performance = function(x,omega){

p = ncol(omega)

A.true = diag(p)

A.x = diag(p)

for (i in 1:p){A.true[i,] = ifelse(abs(omega[i,])<=0,0,1)}

for (i in 1:p){A.x[i,] = ifelse(abs(x[i,])<=0,0,1)}

Edges_I = combn(1:p,2) # Inactive set of ordered pair (i,j)
Edges_x = combn(1:p,2) # Inactive set of ordered pair (i,j)

# Set of Edges

for (t in 1:ncol(Edges_I)){
i = Edges_I[1,t]
j = Edges_I[2,t]
if(A.true[i,j]==0){Edges_I[,t]=c(0,0)}
}

for (t in 1:ncol(Edges_x)){
i = Edges_x[1,t]
j = Edges_x[2,t]
if(A.x[i,j]==0){Edges_x[,t]=c(0,0)}
}

P = which(Edges_I[1,]>0)
N = which(Edges_I[1,]==0)

P.x = which(Edges_x[1,]>0)
N.x = which(Edges_x[1,]==0) 

TP = length(which(P.x %in% P))
TN = length(which(N.x %in% N))

FP = length(P.x) - TP

FN = length(N.x) - TN

TPR = TP/length(P)
FPR = 1 - TN/length(N)

MCC = ((TP*TN)-(FP*FN))/((TP+FP)^0.5*(TP+FN)^0.5*(TN+FP)^0.5*(TN+FN)^0.5)

if (is.nan(MCC)==TRUE){MCC=0}

return(list(TPR=TPR,FPR=FPR,TP=TP,FP=FP))
}

##############################################################
##############################################################
##############################################################

############################################################
##### Compute Adjacency Matrix
############################################################

Adj_mat = function(Omega){

p = dim(Omega)[1]

Adj_mat = matrix(0,p,p)

 for (i in 1:p){
    for (j in 1:p){
      if (abs(Omega[i,j])<=0){Adj_mat[i,j]=0}
      else (Adj_mat[i,j]=1)
    } 
 }
return(Adj_mat)}


############################################################
##### Compute Partial Correlation Matrix
############################################################

# Function that estimates the Regression Weights and Omega

Omega.GGM = function(x,Adj.mat){

x = scale(x)

n = dim(x)[1]

p = dim(x)[2]

Edges_I = combn(1:p,2) # Inactive set of ordered pair (i,j)

# Set of Edges

for (t in 1:ncol(Edges_I)){
i = Edges_I[1,t]
j = Edges_I[2,t]
if(Adj.mat[i,j]==0){Edges_I[,t]=c(0,0)}
}

####  

# Compute Prediction Errors and betas

beta = matrix(0,p,p)

for (i in 1:p){

n_i = c(Edges_I[1,which(Edges_I[2,]==i)],
Edges_I[2,which(Edges_I[1,]==i)])

if (length(n_i)>0){beta[n_i,i] = coef(lm(x[,i] ~ 0 + x[,n_i]))}
}

vareps = x - x%*%beta

# Compute the precision matrix

Omega = matrix(0,p,p)

diag(Omega) = apply(vareps,2,var)^(-1)

for (e in which(Edges_I[2,]>0)){
  i = Edges_I[1,e] 
  j = Edges_I[2,e]
  
  Omega[i,j] = cov(vareps[,i],vareps[,j])*Omega[i,i]*Omega[j,j]
  Omega[j,i] = Omega[i,j] 
}

# Make Omega Positive Definite

lambda_min = eigen(Omega)$values[p]
if (lambda_min < 1e-6){Omega = Omega+(0.1+abs(lambda_min))*diag(p)}

######

R = -cov2cor(Omega)
diag(R) = 1

# Make R Positive Definite

lambda_min = eigen(R)$values[p]
if (lambda_min < 1e-6){R = R+(0.1+abs(lambda_min))*diag(p)}

R = cov2cor(R)

return(R)}

##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################

############################################################
##### True Positives and False Positives
############################################################

# Load data

load(file="sample_1_p_20_n_100.RData")

# Load population partial correlation matrix

load(file="Pcor_1_p_20.RData")

######################################################
######################################################

# Glasso-CV1

load(file = "glasso_cv_1_p_20_n_100.RData")
performance(list.glasso.cv,R)

# PCS-Glasso-CV1

load(file = "pcs_glasso_cv_1_p_20_n_100.RData")
performance(list.glasso.cv.pcs[[2]],R)

######################################################
######################################################

# Glasso-CV1-1se

load(file = "glasso_cv_ls1_1_p_20_n_100.RData")
performance(list.glasso.cv.ls1,R)

# PCS-Glasso-CV1-1se

load(file = "pcs_glasso_cv_ls1_1_p_20_n_100.RData")
performance(list.glasso.cv.ls1.pcs[[2]],R)

######################################################
######################################################

# Glasso-CV2

load(file = "glasso_cv2_1_p_20_n_100.RData")
performance(list.glasso.cv2,R)

# PCS-Glasso-CV2

load(file = "pcs_glasso_cv2_1_p_20_n_100.RData")
performance(list.glasso.cv2.pcs[[2]],R)

######################################################
######################################################

# Glasso-CV2-1se

load(file = "glasso_cv2_ls1_1_p_20_n_100.RData")
performance(list.glasso.cv2.ls1,R)

# PCS-Glasso-CV2-1se

load(file = "pcs_glasso_cv2_ls1_1_p_20_n_100.RData")
performance(list.glasso.cv2.ls1.pcs[[2]],R)

######################################################
######################################################

# Glasso-EBIC

load(file = "glasso_ebic_1_p_20_n_100.RData")
performance(list.glasso.ebic,R)

# PCS-Glasso-EBIC

load(file = "pcs_glasso_ebic_1_p_20_n_100.RData")
performance(list.glasso.ebic.pcs[[2]],R)

######################################################
######################################################

# Glasso-BIC

load(file = "glasso_bic_1_p_20_n_100.RData")
performance(list.glasso.bic,R)

# PCS-Glasso-BIC

load(file = "pcs_glasso_bic_1_p_20_n_100.RData")
performance(list.glasso.bic.pcs[[2]],R)

######################################################
######################################################

# SPACE-CV

load(file = "space_cv_1_p_20_n_100.RData")
performance(list.space.cv[[2]],R)

# PCS-SPACE-CV

load(file = "pcs_space_cv_1_p_20_n_100.RData")
performance(list.space.cv.pcs[[2]],R)

######################################################
######################################################

# SPACE-CV-1se

load(file = "space_cv_ls1_1_p_20_n_100.RData")
performance(list.space.cv.ls1[[2]],R)

# PCS-SPACE-CV-1se

load(file = "pcs_space_cv_ls1_1_p_20_n_100.RData")
performance(list.space.cv.ls1.pcs[[2]],R)

######################################################
######################################################

# SPACE-FSR

load(file = "space_alpha_1_p_20_n_100.RData")
performance(list.space.alpha[[2]],R)

# PCS-SPACE-FSR

load(file = "pcs_space_alpha_1_p_20_n_100.RData")
performance(list.space.alpha.pcs[[2]],R)

######################################################
######################################################

# SPACE-BIC

load(file = "space_bic_1_p_20_n_100.RData")
performance(list.space.bic[[2]],R)

# PCS-SPACE-BIC

load(file = "pcs_space_bic_1_p_20_n_100.RData")
performance(list.space.bic.pcs[[2]],R)

######################################################
######################################################

# NR-AND-CV

load(file = "nei_and_cv_1_p_20_n_100.RData")
performance(list.nei.cv[[2]],R)

# PCS-NR-AND-CV

load(file = "pcs_nei_and_cv_1_p_20_n_100.RData")
performance(list.nei.cv.pcs[[2]],R)

######################################################
######################################################

# NR-AND-CV-1se

load(file = "nei_and_cv_ls1_1_p_20_n_100.RData")
performance(list.nei.cv.ls1[[2]],R)

# PCS-NR-AND-CV-1se

load(file = "pcs_nei_and_cv_ls1_1_p_20_n_100.RData")
performance(list.nei.cv.ls1.pcs[[2]],R)

######################################################
######################################################

# NR-AND-FSR

load(file = "nei_and_alpha_1_p_20_n_100.RData")
performance(list.nei.alpha[[2]],R)

# PCS-NR-AND-FSR

load(file = "pcs_nei_and_alpha_1_p_20_n_100.RData")
performance(list.nei.alpha.pcs[[2]],R)

######################################################
######################################################

# NR-AND-BIC

load(file = "nei_and_bic_1_p_20_n_100.RData")
performance(list.nei.bic[[2]],R)

# PCS-NR-AND-BIC

load(file = "pcs_nei_and_bic_1_p_20_n_100.RData")
performance(list.nei.bic.pcs[[2]],R)

######################################################
######################################################

# NR-OR-CV

load(file = "nei_or_cv_1_p_20_n_100.RData")
performance(list.nei.cv[[2]],R)

# PCS-NR-OR-CV

load(file = "pcs_nei_or_cv_1_p_20_n_100.RData")
performance(list.nei.cv.pcs[[2]],R)

######################################################
######################################################

# NR-OR-CV-1se

load(file = "nei_or_cv_ls1_1_p_20_n_100.RData")
performance(list.nei.cv.ls1[[2]],R)

# PCS-NR-OR-CV-1se

load(file = "pcs_nei_or_cv_ls1_1_p_20_n_100.RData")
performance(list.nei.cv.ls1.pcs[[2]],R)

######################################################
######################################################

# NR-OR-FSR

load(file = "nei_or_alpha_1_p_20_n_100.RData")
performance(list.nei.alpha[[2]],R)

# PCS-NR-OR-FSR

load(file = "pcs_nei_or_alpha_1_p_20_n_100.RData")
performance(list.nei.alpha.pcs[[2]],R)

######################################################
######################################################

# NR-OR-BIC

load(file = "nei_or_bic_1_p_20_n_100.RData")
performance(list.nei.bic[[2]],R)

# PCS-NR-OR-BIC

load(file = "pcs_nei_or_bic_1_p_20_n_100.RData")
performance(list.nei.bic.pcs[[2]],R)

######################################################
######################################################

# Ridge-CV

load(file = "ridge_cv_1_p_20_n_100.RData")
performance(list.ridge.cv[[2]],R)

# PCS-Ridge-CV

load(file = "pcs_ridge_cv_1_p_20_n_100.RData")
performance(list.ridge.cv.pcs[[2]],R)

######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################

############################################################
##### HeatMaps
############################################################

par(mfrow=c(4,5))

corrplot(Adj_mat(R),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "True", mar=c(0,0,1,0))

# Glasso-CV1

load(file = "glasso_cv_1_p_20_n_100.RData")
p.glasso.cv1 = corrplot(Adj_mat(list.glasso.cv),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "Glasso-CV1", mar=c(0,0,1,0))

# Glasso-CV1-1se

load(file = "glasso_cv_ls1_1_p_20_n_100.RData")
p.glasso.cv1.1se = corrplot(Adj_mat(list.glasso.cv.ls1),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "Glasso-CV1-1se", mar=c(0,0,1,0))

# Glasso-CV2

load(file = "glasso_cv2_1_p_20_n_100.RData")
p.glasso.cv2 = corrplot(Adj_mat(list.glasso.cv2),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "Glasso-CV2", mar=c(0,0,1,0))

# Glasso-CV2-1se

load(file = "glasso_cv2_ls1_1_p_20_n_100.RData")
p.glasso.cv2.1se = corrplot(Adj_mat(list.glasso.cv2.ls1),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "Glasso-CV2-1se", mar=c(0,0,1,0))

# Glasso-EBIC

load(file = "glasso_ebic_1_p_20_n_100.RData")
p.glasso.ebic = corrplot(Adj_mat(list.glasso.ebic),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "Glasso-EBIC", mar=c(0,0,1,0))

# Glasso-BIC

load(file = "glasso_bic_1_p_20_n_100.RData")
p.glasso.bic = corrplot(Adj_mat(list.glasso.bic),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "Glasso-BIC", mar=c(0,0,1,0))

# SPACE-CV

load(file = "space_cv_1_p_20_n_100.RData")
p.space.cv = corrplot(Adj_mat(list.space.cv[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "SPACE-CV", mar=c(0,0,1,0))

# SPACE-CV-1se

load(file = "space_cv_ls1_1_p_20_n_100.RData")
p.space.cv.1se = corrplot(Adj_mat(list.space.cv.ls1[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "SPACE-CV-1se", mar=c(0,0,1,0))

# SPACE-FSR

load(file = "space_alpha_1_p_20_n_100.RData")
p.space.fsr = corrplot(Adj_mat(list.space.alpha[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "SPACE-FSR", mar=c(0,0,1,0))

# SPACE-BIC

load(file = "space_bic_1_p_20_n_100.RData")
p.space.bic = corrplot(Adj_mat(list.space.bic[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "SPACE-BIC", mar=c(0,0,1,0))

# NR-AND-CV

load(file = "nei_and_cv_1_p_20_n_100.RData")
p.nr.and.cv = corrplot(Adj_mat(list.nei.cv[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "NR-AND-CV", mar=c(0,0,1,0))

# NR-AND-CV-1se

load(file = "nei_and_cv_ls1_1_p_20_n_100.RData")
p.nr.and.cv.1se = corrplot(Adj_mat(list.nei.cv.ls1[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "NR-AND-CV-1se", mar=c(0,0,1,0))

# NR-AND-FSR

load(file = "nei_and_alpha_1_p_20_n_100.RData")
p.nr.and.fsr = corrplot(Adj_mat(list.nei.alpha[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "NR-AND-FSR", mar=c(0,0,1,0))

# NR-AND-BIC

load(file = "nei_and_bic_1_p_20_n_100.RData")
p.nr.and.bic = corrplot(Adj_mat(list.nei.bic[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "NR-AND-BIC", mar=c(0,0,1,0))

# NR-OR-CV

load(file = "nei_or_cv_1_p_20_n_100.RData")
p.nr.or.cv = corrplot(Adj_mat(list.nei.cv[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "NR-OR-CV", mar=c(0,0,1,0))

# NR-OR-CV-1se

load(file = "nei_or_cv_ls1_1_p_20_n_100.RData")
p.nr.or.cv.1se = corrplot(Adj_mat(list.nei.cv.ls1[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "NR-OR-CV-1se", mar=c(0,0,1,0))

# NR-OR-FSR

load(file = "nei_or_alpha_1_p_20_n_100.RData")
p.nr.or.fsr = corrplot(Adj_mat(list.nei.alpha[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "NR-OR-FSR", mar=c(0,0,1,0))

# NR-OR-BIC

load(file = "nei_or_bic_1_p_20_n_100.RData")
p.nr.or.bic = corrplot(Adj_mat(list.nei.bic[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "NR-OR-BIC", mar=c(0,0,1,0))

# Ridge-CV

load(file = "ridge_cv_1_p_20_n_100.RData")
p.ridge.cv = corrplot(Adj_mat(list.ridge.cv[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),,
title = "Ridge-CV", mar=c(0,0,1,0))

######################################################
######################################################
######################################################
######################################################
######################################################
######################################################

par(mfrow=c(4,5))

corrplot(Adj_mat(R),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "True", mar=c(0,0,1,0))

# PCS-Glasso-CV1

load(file = "pcs_glasso_cv_1_p_20_n_100.RData")
p.glasso.cv1.pcs = corrplot(Adj_mat(list.glasso.cv.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-Glasso-CV1", mar=c(0,0,1,0))

# PCS-Glasso-CV1-1se

load(file = "pcs_glasso_cv_ls1_1_p_20_n_100.RData")
p.glasso.cv1.1se.pcs = corrplot(Adj_mat(list.glasso.cv.ls1.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-Glasso-CV1-1se", mar=c(0,0,1,0))

# PCS-Glasso-CV2

load(file = "pcs_glasso_cv2_1_p_20_n_100.RData")
p.glasso.cv2.pcs = corrplot(Adj_mat(list.glasso.cv2.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-Glasso-CV2", mar=c(0,0,1,0))

# PCS-Glasso-CV2-1se

load(file = "pcs_glasso_cv2_ls1_1_p_20_n_100.RData")
p.glasso.cv2.1se.pcs = corrplot(Adj_mat(list.glasso.cv2.ls1.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-Glasso-CV2-1se", mar=c(0,0,1,0))

# PCS-Glasso-EBIC

load(file = "pcs_glasso_ebic_1_p_20_n_100.RData")
p.glasso.ebic.pcs = corrplot(Adj_mat(list.glasso.ebic.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-Glasso-EBIC", mar=c(0,0,1,0))

# PCS-Glasso-BIC

load(file = "pcs_glasso_bic_1_p_20_n_100.RData")
p.glasso.bic.pcs = corrplot(Adj_mat(list.glasso.bic.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-Glasso-BIC", mar=c(0,0,1,0))

# PCS-SPACE-CV

load(file = "pcs_space_cv_1_p_20_n_100.RData")
p.space.cv.pcs = corrplot(Adj_mat(list.space.cv.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-SPACE-CV", mar=c(0,0,1,0))

# PCS-SPACE-CV-1se

load(file = "pcs_space_cv_ls1_1_p_20_n_100.RData")
p.space.cv.1se.pcs = corrplot(Adj_mat(list.space.cv.ls1.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-SPACE-CV-1se", mar=c(0,0,1,0))

# PCS-SPACE-FSR

load(file = "pcs_space_alpha_1_p_20_n_100.RData")
p.space.fsr.pcs = corrplot(Adj_mat(list.space.alpha.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-SPACE-FSR", mar=c(0,0,1,0))

# PCS-SPACE-BIC

load(file = "pcs_space_bic_1_p_20_n_100.RData")
p.space.bic.pcs = corrplot(Adj_mat(list.space.bic.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-SPACE-BIC", mar=c(0,0,1,0))

# PCS-NR-AND-CV

load(file = "pcs_nei_and_cv_1_p_20_n_100.RData")
p.nr.and.cv.pcs = corrplot(Adj_mat(list.nei.cv.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-NR-AND-CV", mar=c(0,0,1,0))

# PCS-NR-AND-CV-1se

load(file = "pcs_nei_and_cv_ls1_1_p_20_n_100.RData")
p.nr.and.cv.1se.pcs = corrplot(Adj_mat(list.nei.cv.ls1.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-NR-AND-CV-1se", mar=c(0,0,1,0))

# PCS-NR-AND-FSR

load(file = "pcs_nei_and_alpha_1_p_20_n_100.RData")
p.nr.and.fsr.pcs = corrplot(Adj_mat(list.nei.alpha.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-NR-AND-FSR", mar=c(0,0,1,0))

# PCS-NR-AND-BIC

load(file = "pcs_nei_and_bic_1_p_20_n_100.RData")
p.nr.and.bic.pcs = corrplot(Adj_mat(list.nei.bic.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-NR-AND-BIC", mar=c(0,0,1,0))

# PCS-NR-OR-CV

load(file = "pcs_nei_or_cv_1_p_20_n_100.RData")
p.nr.or.cv.pcs = corrplot(Adj_mat(list.nei.cv.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-NR-OR-CV", mar=c(0,0,1,0))

# PCS-NR-OR-CV-1se

load(file = "pcs_nei_or_cv_ls1_1_p_20_n_100.RData")
p.nr.or.cv.1se.pcs = corrplot(Adj_mat(list.nei.cv.ls1.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-NR-OR-CV-1se", mar=c(0,0,1,0))

# PCS-NR-OR-FSR

load(file = "pcs_nei_or_alpha_1_p_20_n_100.RData")
p.nr.or.fsr.pcs = corrplot(Adj_mat(list.nei.alpha.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-NR-OR-FSR", mar=c(0,0,1,0))

# PCS-NR-OR-BIC

load(file = "pcs_nei_or_bic_1_p_20_n_100.RData")
p.nr.or.bic.pcs = corrplot(Adj_mat(list.nei.bic.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-NR-OR-BIC", mar=c(0,0,1,0))

# PCS-Ridge-CV

load(file = "pcs_ridge_cv_1_p_20_n_100.RData")
p.ridge.cv.pcs = corrplot(Adj_mat(list.ridge.cv.pcs[[2]]),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200),
title = "PCS-Ridge-CV", mar=c(0,0,1,0))

######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################

############################################################
##### Networks
############################################################

par(mfrow=c(4,5))

qgraph(Omega.GGM(x,Adj_mat(R)), graph = "cor", fade=FALSE)
title("True",line=3)

# Glasso-CV1

load(file = "glasso_cv_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.glasso.cv)), graph = "cor", fade=FALSE)
title("Glasso-CV1",line=3)

# Glasso-CV1-1se

load(file = "glasso_cv_ls1_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.glasso.cv.ls1)), graph = "cor", fade=FALSE)
title("Glasso-CV1-1se",line=3)

# Glasso-CV2

load(file = "glasso_cv2_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.glasso.cv2)), graph = "cor", fade=FALSE)
title("Glasso-CV2",line=3)

# Glasso-CV2-1se

load(file = "glasso_cv2_ls1_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.glasso.cv2.ls1)), graph = "cor", fade=FALSE)
title("Glasso-CV2-1se",line=3)

# Glasso-EBIC

load(file = "glasso_ebic_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.glasso.ebic)), graph = "cor", fade=FALSE)
title("Glasso-EBIC",line=3)

# Glasso-BIC

load(file = "glasso_bic_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.glasso.bic)), graph = "cor", fade=FALSE)
title("Glasso-BIC",line=3)

# SPACE-CV

load(file = "space_cv_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.space.cv[[2]])), graph = "cor", fade=FALSE)
title("SPACE-CV",line=3)

# SPACE-CV-1se

load(file = "space_cv_ls1_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.space.cv.ls1[[2]])), graph = "cor", fade=FALSE)
title("SPACE-CV-1se",line=3)

# SPACE-FSR

load(file = "space_alpha_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.space.alpha[[2]])), graph = "cor", fade=FALSE)
title("SPACE-FSR",line=3)

# SPACE-BIC

load(file = "space_bic_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.space.bic[[2]])), graph = "cor", fade=FALSE)
title("SPACE-BIC",line=3)

# NR-AND-CV

load(file = "nei_and_cv_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.cv[[2]])), graph = "cor", fade=FALSE)
title("NR-AND-CV",line=3)

# NR-AND-CV-1se

load(file = "nei_and_cv_ls1_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.cv.ls1[[2]])), graph = "cor", fade=FALSE)
title("NR-AND-CV-1se",line=3)

# NR-AND-FSR

load(file = "nei_and_alpha_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.alpha[[2]])), graph = "cor", fade=FALSE)
title("NR-AND-FSR",line=3)

# NR-AND-BIC

load(file = "nei_and_bic_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.bic[[2]])), graph = "cor", fade=FALSE)
title("NR-AND-BIC",line=3)

# NR-OR-CV

load(file = "nei_or_cv_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.cv[[2]])), graph = "cor", fade=FALSE)
title("NR-OR-CV",line=3)

# NR-OR-CV-1se

load(file = "nei_or_cv_ls1_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.cv.ls1[[2]])), graph = "cor", fade=FALSE)
title("NR-OR-CV-1se",line=3)

# NR-OR-FSR

load(file = "nei_or_alpha_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.alpha[[2]])), graph = "cor", fade=FALSE)
title("NR-OR-FSR",line=3)

# NR-OR-BIC

load(file = "nei_or_bic_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.bic[[2]])), graph = "cor", fade=FALSE)
title("NR-OR-BIC",line=3)

# Ridge-CV

load(file = "ridge_cv_1_p_20_n_100.RData")
qgraph(solve(list.ridge.cv[[2]]), graph = "pcor", fade=FALSE)
title("Ridge-CV",line=3)

######################################################
######################################################
######################################################
######################################################

par(mfrow=c(4,5))

qgraph(Omega.GGM(x,Adj_mat(R)), graph = "cor", fade=FALSE)
title("True",line=3)

# PCS-Glasso-CV1

load(file = "pcs_glasso_cv_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.glasso.cv.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-Glasso-CV1",line=3)

# PCS-Glasso-CV1-1se

load(file = "pcs_glasso_cv_ls1_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.glasso.cv.ls1.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-Glasso-CV1-1se",line=3)

# PCS-Glasso-CV2

load(file = "pcs_glasso_cv2_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.glasso.cv2.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-Glasso-CV2",line=3)

# PCS-Glasso-CV2-1se

load(file = "pcs_glasso_cv2_ls1_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.glasso.cv2.ls1.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-Glasso-CV2-1se",line=3)

# PCS-Glasso-EBIC

load(file = "pcs_glasso_ebic_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.glasso.ebic.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-Glasso-EBIC",line=3)

# PCS-Glasso-BIC

load(file = "pcs_glasso_bic_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.glasso.bic.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-Glasso-BIC",line=3)

# PCS-SPACE-CV

load(file = "pcs_space_cv_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.space.cv.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-SPACE-CV",line=3)

# PCS-SPACE-CV-1se

load(file = "pcs_space_cv_ls1_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.space.cv.ls1.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-SPACE-CV-1se",line=3)

# PCS-SPACE-FSR

load(file = "pcs_space_alpha_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.space.alpha.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-SPACE-FSR",line=3)

# PCS-SPACE-BIC

load(file = "pcs_space_bic_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.space.bic.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-SPACE-BIC",line=3)

# PCS-NR-AND-CV

load(file = "pcs_nei_and_cv_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.cv.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-NR-AND-CV",line=3)

# PCS-NR-AND-CV-1se

load(file = "pcs_nei_and_cv_ls1_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.cv.ls1.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-NR-AND-CV-1se",line=3)

# PCS-NR-AND-FSR

load(file = "pcs_nei_and_alpha_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.alpha.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-NR-AND-FSR",line=3)

# PCS-NR-AND-BIC

load(file = "pcs_nei_and_bic_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.bic.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-NR-AND-BIC",line=3)

# PCS-NR-OR-CV

load(file = "pcs_nei_or_cv_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.cv.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-NR-OR-CV",line=3)

# PCS-NR-OR-CV-1se

load(file = "pcs_nei_or_cv_ls1_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.cv.ls1.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-NR-OR-CV-1se",line=3)

# PCS-NR-OR-FSR

load(file = "pcs_nei_or_alpha_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.alpha.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-NR-OR-FSR",line=3)

# PCS-NR-OR-BIC

load(file = "pcs_nei_or_bic_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.nei.bic.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-NR-OR-BIC",line=3)

# PCS-Ridge-CV

load(file = "pcs_ridge_cv_1_p_20_n_100.RData")
qgraph(Omega.GGM(x,Adj_mat(list.ridge.cv.pcs[[2]])), graph = "cor", fade=FALSE)
title("PCS-Ridge-CV",line=3)


##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################

############################################################
##### PCS-Networks Combination 
############################################################

# Create a matrix that combines the 19 PCS estimates

p = ncol(x)

Adj_mat_average = matrix(0,p,p)

#############################################################

load(file = "PCS_glasso_cv_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat(list.glasso.cv.pcs[[2]])

#############################################################

load(file = "PCS_glasso_cv_ls1_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.glasso.cv.ls1.pcs[[2]])

#############################################################

load(file = "PCS_glasso_cv2_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat(list.glasso.cv.pcs[[2]])

#############################################################

load(file = "PCS_glasso_cv2_ls1_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.glasso.cv2.ls1.pcs[[2]])

#############################################################

load(file = "PCS_glasso_bic_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.glasso.bic.pcs[[2]])

#############################################################

load(file = "PCS_glasso_ebic_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.glasso.ebic.pcs[[2]])

#############################################################

load(file = "PCS_space_cv_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.space.cv.pcs[[2]])

#############################################################

load(file = "PCS_space_cv_ls1_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.space.cv.ls1.pcs[[2]])

#############################################################

load(file = "PCS_space_alpha_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.space.alpha.pcs[[2]])

#############################################################

load(file = "PCS_space_bic_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.space.bic.pcs[[2]])

#############################################################

load(file = "PCS_nei_and_cv_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.nei.cv.pcs[[2]])

#############################################################

load(file = "PCS_nei_and_cv_ls1_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.nei.cv.ls1.pcs[[2]])

#############################################################

load(file = "PCS_nei_and_alpha_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.nei.alpha.pcs[[2]])

#############################################################

load(file = "PCS_nei_and_bic_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.nei.bic.pcs[[2]])

#############################################################

load(file = "PCS_nei_or_cv_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.nei.cv.pcs[[2]])

#############################################################

load(file = "PCS_nei_or_cv_ls1_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.nei.cv.ls1.pcs[[2]])

#############################################################

load(file = "PCS_nei_or_alpha_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.nei.alpha.pcs[[2]])

#############################################################

load(file = "PCS_nei_or_bic_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.nei.bic.pcs[[2]])

#############################################################

load(file = "PCS_ridge_cv_1_p_20_n_100.RData")
Adj_mat_average = Adj_mat_average + Adj_mat(list.ridge.cv.pcs[[2]])

#############################################################
#############################################################
#############################################################
#############################################################

# Generate a matrix that contains the edges with a frequency of appearance of at least in two of the estimation procedures

Adj_mat_average = Adj_mat_average

for (i in 1:p){
for (j in 1:p){
if (abs(Adj_mat_average[i,j])<=1){Adj_mat_average[i,j]=0}
else (Adj_mat_average[i,j]=1)
}}

# Compute Partial Correlation MAtrix with the combination of the models

R.combination = Omega.GGM(x,Adj_mat_average)

# Network plot

qgraph(R.combination, graph = "cor", layout = "circular", fade=FALSE)

# Heatmap

corrplot(Adj_mat(R.combination),method='color',is.corr = FALSE,
cl.lim=c(0,1), col=colorRampPalette(c("blue","white","black"))(200))

