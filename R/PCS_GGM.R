#################################################################
##########         R Code by Ginette Lafit        ###############
##########         ginette.lafit@kuleuven.be      ###############
#################################################################

#################################################################
#################################################################
#################################################################
##### GGM Estimation using Partial Correlation Screening (PCS)
#################################################################
#################################################################
#################################################################

# Function that estimates the Regression Weights and Omega using Partial Correlation Screening (PCS)
# the input is the (n x p) data matrix (i.e., x), the estimated partial correlation matrix and the threshold tau
# n is the number of rows of x
# p is the number of columns of x
# R is the estimated partial correlation matrix 
# tau is the threshold to be used by the PCS approach. tau is included in the interval [0,1]

# The function returns the estimated precision matrix and a matrix that contains the regression weights for each node 
# in the graph, as a function of the elements in the neighborhood

PCS_GGM = function(x,R,tau){

x = scale(x)

n = dim(x)[1]

p = dim(x)[2]

Edges_I = combn(1:p,2) # Inactive set of ordered pair (i,j)

# Apply the PCS approach to estimate the set of edges by thresholding the elements of the matrix R

for (i in 1:p){R[i,] = ifelse(abs(R[i,])<=tau,0,R[i,])}

for (t in 1:ncol(Edges_I)){
i = Edges_I[1,t]
j = Edges_I[2,t]
if(R[i,j]==0){Edges_I[,t]=c(0,0)}
}

####  

# Compute prediction errors and regression weights 

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

# Make Omega Positive Definite by adding to the diagonal a constant equal to 0.1 + abs(lambda_min)
# where lambda_min is the minimum eigenvalue of Omega

lambda_min = eigen(Omega)$values[p]
if (lambda_min < 1e-6){Omega = Omega+(0.1+abs(lambda_min))*diag(p)}

######
return(list(beta,Omega))}