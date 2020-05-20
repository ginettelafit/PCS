#################################################################
##########         R Code by Ginette Lafit        ###############
##########         ginette.lafit@kuleuven.be      ###############
#################################################################

###########################################################
###########################################################
###########################################################
##### Nodewise regression BIC
###########################################################
###########################################################
###########################################################

# Nodewise regression + BIC

# Function that estimates the regularization parameter in nodewise regression by minimizing BIC
# the input is the (n x p) data matrix (i.e., x)
# n is the number of rows of x
# p is the number of columns of x

# The function returns the regularization parameter that minimizes BIC 

BIC_Neigh = function(x){

library(glmnet)

x = scale(x)
 
n = nrow(x)
p = ncol(x)

lambda.opt = rep(0,p)

for (i in 1:p){
fit.neigh = glmnet(x[,-i], x[,i], family='gaussian',nlambda=100)
beta.list = as.matrix(fit.neigh$beta)
lambdalist = fit.neigh$lambda
BIC = rep(0,ncol(beta.list))
for (k in 1:ncol(beta.list)){
BIC[k] = n*log(sum((x[,i] - x[,-i]%*%beta.list[,k])^2)) + log(n)*sum(abs(beta.list[,k])>0)
}
ind = which.min(BIC)
lambda.opt[i] = lambdalist[ind]
}

return(lambda.opt)  
}