#################################################################
##########         R Code by Ginette Lafit        ###############
##########         ginette.lafit@kuleuven.be      ###############
#################################################################

###########################################################
###########################################################
###########################################################
##### SPACE BIC
###########################################################
###########################################################
###########################################################

# SPACE + BIC

# Function that selects the regularization parameter on SPACE that minimizes BIC

# the input is the (n x p) data matrix (i.e., x)
# n is the number of rows of x
# p is the number of columns of x

# The function returns the regularization parameter that minimizes BIC

BIC_space = function(x){

library(space)

x = scale(x)
 
n = nrow(x)
p = ncol(x)

l_max = sqrt(n)*qnorm(1-0.0001/(2*p^2))

l_min = sqrt(n)*qnorm(1-0.9/(2*p^2))

lamlist = seq(l_min,l_max,length=100)

loss.re = rep(0,length(lamlist))

vareps = lapply(1:length(lamlist),function(i) 
space.joint(x, lam1=lamlist[[i]], lam2=0, weight=NULL, iter=3)[[1]])

eps = lapply(1:length(lamlist),function(i) PCS_GGM(x,vareps[[i]],0)[[1]])

for (k in 1:length(lamlist)){
beta = eps[[k]]
BIC = rep(0,p)
for (i in 1:p){
BIC[i] = n*log(sum((x[,i] - x[,-i]%*%beta[i,-i])^2)) + log(n)*sum(abs(beta[i,-i])>0)
}
loss.re[k] = sum(BIC)
} 

ind = which.min(loss.re)
lamlist_opt = lamlist[ind]

res = lamlist_opt
return(res)  
}

