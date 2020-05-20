#################################################################
##########         R Code by Ginette Lafit        ###############
##########         ginette.lafit@kuleuven.be      ###############
#################################################################

###########################################################
###########################################################
###########################################################
##### Nodewise regression k-fold cross-validation
###########################################################
###########################################################
###########################################################

# Nodewise regression + CV

# Function that estimates the regularization parameter in nodewise regression using k-fold cross-validation
# the input is the (n x p) data matrix (i.e., x)
# n is the number of rows of x
# p is the number of columns of x
# fold is the number of folds

# The function returns 2 regularization parameters: lambda.opt and lambda.ls1
# lambda.opt minimizes the mean squared error
# lambda.ls1 is the regularization parameter using the one-standard-error-rule

cv_Neigh = function(x, fold){

library(glmnet)

x = scale(x)

p = ncol(x)

lambda.fit = lapply(1:p, function(i) cv.glmnet(x[,-i], x[,i], family='gaussian', nfolds=fold))

lambda.opt = unlist(lambda.fit)

lambda.opt = rep(0,p)
lambda.lse = rep(0,p)

for (i in 1:p){
lambda.opt[i] = lambda.fit[[i]]$lambda.min
lambda.lse[i] = lambda.fit[[i]]$lambda.1se
}

######

return(list(lambda.opt,lambda.lse))}