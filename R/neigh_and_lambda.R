#################################################################
##########         R Code by Ginette Lafit        ###############
##########         ginette.lafit@kuleuven.be      ###############
#################################################################

###########################################################
###########################################################
###########################################################
##### GGM estimation using Nodewise Regression 'AND' Rule
###########################################################
###########################################################
###########################################################

# Estimate the GGM using Nodewise regression + AND rule

# Function that estimates the Regression Weights and the Precision Matrix using nodewise regressions and the AND rule
# the input is the (n x p) data matrix (i.e., x) and the regularization parameter lambda
# n is the number of rows of x
# p is the number of columns of x
# lambda is the regularization parameter use to estimate the nodewise regressions using lasso (i.e., lambda>0)

# The function returns a matrix beta that contains the regression weights for each node 
# in the graph, as a function of the elements in the neighborhood
# and the precision matrix Omega.hat

neigh_and_lambda = function(x,lambda){

library(glmnet)

x = scale(x)

p = ncol(x)

n = nrow(x)

Omega = matrix(0,p,p)

for (i in 1:p){
Omega[i,-i] = as.matrix(glmnet(x[,-i], x[,i], family='gaussian', lambda = lambda[i])$beta)
}

for (i in 1:p){
for (j in 1:p){
Omega[i,j] = min(abs(Omega[i,j]),abs(Omega[j,i]))
Omega[j,i] = Omega[i,j] 
}}

# Compute the precision matrix

Edges_I = combn(1:p,2) # Inactive set of ordered pair (i,j)

# Candidate Set of Edges

for (t in 1:ncol(Edges_I)){
i = Edges_I[1,t]
j = Edges_I[2,t]
if(Omega[i,j]==0){Edges_I[,t]=c(0,0)}
}

# Compute betas

beta = matrix(0,p,p)

for (i in 1:p){

n_i = c(Edges_I[1,which(Edges_I[2,]==i)],
Edges_I[2,which(Edges_I[1,]==i)])

if (length(n_i)>0){beta[n_i,i] = coef(lm(x[,i] ~ 0 + x[,n_i]))}
}

vareps = x - x%*%beta

# Estimate Precision Matrix

Omega.hat = matrix(0,p,p)

diag(Omega.hat) = apply(vareps,2,var)^(-1)

for (e in which(Edges_I[2,]>0)){
  i = Edges_I[1,e] 
  j = Edges_I[2,e]
  Omega.hat[i,j] = cov(vareps[,i],vareps[,j])*Omega.hat[i,i]*Omega.hat[j,j]
  Omega.hat[j,i] = Omega.hat[i,j] 
}

# Make Omega Positive Definite by adding to the diagonal a constant equal to 0.1 + abs(lambda_min)
# where lambda_min is the minimum eigenvalue of Omega

lambda_min = eigen(Omega.hat)$values[p]
if (lambda_min < 1e-6){Omega.hat = Omega.hat+(0.1+abs(lambda_min))*diag(p)}

######

return(list(beta,Omega.hat))}
