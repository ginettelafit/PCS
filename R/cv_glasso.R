#################################################################
##########         R Code by Ginette Lafit        ###############
##########         ginette.lafit@kuleuven.be      ###############
#################################################################

###########################################################
###########################################################
###########################################################
##### Glasso k-fold cross-validation
###########################################################
###########################################################
###########################################################

# Glasso + CV

# Function that selects the regularization parameter on Glasso by perfoming k-fold cross-validation

# the input is the (n x p) data matrix (i.e., x)
# n is the number of rows of x
# p is the number of columns of x
# fold is the number of folds

# The function returns 4 regularization parameters: alpha_opt and alpha_ls1
# rho_opt minimizes the negative log-likelihood
# rho_ls1 is the regularization parameter using the one-standard-error-rule when the loss function is the negative log-likelihood
# rho2_opt minimizes the mean squared error
# rho2_ls1 is the regularization parameter using the one-standard-error-rule when the loss function is the mean squared error


cv_glasso = function(x,fold){

library(huge)
library(glasso)
 
  cv.part = function(n, k) {
    ntest = floor(n/k)
    ntrain = n - ntest
    ind = sample(n)
    trainMat = matrix(NA, nrow=ntrain, ncol=k)
    testMat = matrix(NA, nrow=ntest, ncol=k)
    nn = 1:n
    for (j in 1:k) {
      sel = ((j-1)*ntest+1):(j*ntest)
      testMat[,j] = ind[sel ]
      sel2 = nn[ !(nn %in% sel) ]
      trainMat[,j] = ind[sel2]
    }
    return(list(trainMat=trainMat,testMat=testMat))
  }
  
  loss_likelihood = function(Sigma, Omega){
    tmp = (sum(diag(Sigma%*%Omega)) - log(det(Omega)) - dim(Omega)[1])
    if(is.finite(tmp)) return(tmp)
    else tmp = Inf
    return(tmp)
  }

  beta_mat = function(Omega){
    p = ncol(Omega)
    beta = matrix(0,p,p)
    for (i in 1:p){beta[,i] = -Omega[,i]/Omega[i,i]}
    diag(beta) = 0
    return(beta)
  }

  loss.cv = function(x.train,x.test,rholist){
    Omegalist = huge(x.train,lambda = rholist ,scr = F, method = "glasso",nlambda=100)$icov
    loss.re = unlist(lapply(1:100, function(i) loss_likelihood(cov(x.test),Omegalist[[i]])))
    loss2.re = unlist(lapply(1:100, function(i) sum(colSums((x.test - x.test%*%beta_mat(Omegalist[[i]]))^2))))
  return(list(loss.re,loss2.re))
  } 
  
  n = nrow(x)
  part.list = cv.part(n, fold)
  
  rholist = seq(0.001,max(abs(cov(x))),length=100)   

  loss.list = lapply(1:fold, function(k) 
  loss.cv(x[part.list$trainMat[, k], ],x[part.list$testMat[, k], ],rholist))

  loss.re = matrix(unlist(lapply(1:fold, function(k) loss.list[[k]][[1]])),ncol = fold, byrow = FALSE)
  loss2.re = matrix(unlist(lapply(1:fold, function(k) loss.list[[k]][[2]])),ncol = fold, byrow = FALSE)
  
  loss.mean = apply(loss.re, 1, mean)
  std.error = apply(loss.re, 1, sd)/sqrt(fold)
  ind = which.min(loss.mean)

  loss.mean.ls = loss.mean[ind:length(loss.mean)]
  rho.list.ls = rholist[ind:length(loss.mean)]

  loss.mean.ls.max = loss.mean.ls[loss.mean.ls <= loss.mean[ind] + std.error[ind]]
  rho.list.ls.max = rho.list.ls[loss.mean.ls <= loss.mean[ind] + std.error[ind]]
  ind.ls = length(loss.mean.ls.max)

  rho_opt = rholist[ind]
  rho_ls1 = rho.list.ls.max[ind.ls]

  ###

  loss2.mean = apply(loss2.re, 1, mean)
  std2.error = apply(loss2.re, 1, sd)/sqrt(fold)
  ind2 = which.min(loss2.mean) 

  loss2.mean.ls = loss2.mean[ind2:length(loss2.mean)]
  rho2.list.ls = rholist[ind2:length(loss2.mean)]

  loss2.mean.ls.max = loss2.mean.ls[loss2.mean.ls <= loss2.mean[ind2] + std2.error[ind2]]
  rho2.list.ls.max = rho2.list.ls[loss2.mean.ls <= loss2.mean[ind2] + std2.error[ind2]]
  ind2.ls = length(loss2.mean.ls.max)

  rho2_opt = rholist[ind2]
  rho2_ls1 = rho2.list.ls.max[ind2.ls]
  
  res = list(rho_opt=rho_opt,rho_ls1=rho_ls1,rho2_opt=rho2_opt,rho2_ls1=rho2_ls1)
  return(res)  
}