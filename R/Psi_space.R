#################################################################
##########         R Code by Ginette Lafit        ###############
##########         ginette.lafit@kuleuven.be      ###############
#################################################################

###########################################################
###########################################################
###########################################################
##### SPACE
###########################################################
###########################################################
###########################################################

# SPACE 

# Function that estimates the Partial Correlation Matrix using SPACE
# the input is the (n x p) data matrix (i.e., x), the regularization parameter lambda.opt
# n is the number of rows of x
# p is the number of columns of x
# lambda.opt is the regularization parameter (i.e., lambda.opt>0)

# The function returns the estimated partial correlation matrix

Psi_space = function(x,lambda.opt){

library(space)

#Initialization

x = scale(x)

n = dim(x)[1]

p = dim(x)[2]

fit = space.joint(x, lam1=lambda.opt, lam2=0, weight=NULL, iter=3)

R = fit[[1]]

return(R)}