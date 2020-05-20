#################################################################
##########         R Code by Ginette Lafit        ###############
##########         ginette.lafit@kuleuven.be      ###############
#################################################################

###########################################################
###########################################################
###########################################################
##### Ridge regression k-fold cross-validation
###########################################################
###########################################################
###########################################################

# Ridge regression + CV

# Function that estimates the regularization parameter in ridge regression using k-fold cross-validation
# the input is the (n x p) data matrix (i.e., x)
# n is the number of rows of x
# p is the number of columns of x
# fold is the number of folds

# The function returns the regularization parameters that minimizes the mean squared error

cv_ridge = function(x, fold){

library(glmnet)

x = scale(x)

p = ncol(x)

lambda.opt = lapply(1:p, function(i) cv.glmnet(x[,-i], x[,i], family='gaussian',alpha = 0, nfolds=fold)$lambda.min)

lambda.opt = unlist(lambda.opt)

######

return(lambda.opt)}