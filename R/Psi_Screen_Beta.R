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

# Function that estimates the Regression Weights using Partial Correlation Screening (PCS)
# the input is the (n x p) data matrix (i.e., x), the estimated partial correlation matrix and the threshold tau
# n is the number of rows of x
# p is the number of columns of x
# R is the estimated partial correlation matrix 
# tau is the threshold to be used by the PCS approach. tau is included in the interval [0,1]

# The function returns a matrix that contains the regression weights for each node 
# in the graph, as a function of the elements in the neighborhood

Psi_Screen_Beta = function(x,R,tau){

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

# Compute Neighborhood

nei_i = function(i,R){
ne_i = which(abs(R[i,])>0)
ne_i = ne_i[-which(ne_i==i)]
return(ne_i)}

Nei = lapply(1:p, function(i) nei_i(i,R))

####  

# Compute Prediction Errors and betas

beta_i = function(x,i,ne_i){
b_i = rep(0,ncol(x))
if (length(ne_i)>0){b_i[ne_i] = coef(lm(x[,i] ~ 0 + x[,ne_i]))}
return(b_i)}

beta_list = lapply(1:p, function(i) beta_i(x,i,Nei[[i]]))

beta = t(matrix(unlist(beta_list), ncol = p, byrow = TRUE))

######

return(list(beta))}