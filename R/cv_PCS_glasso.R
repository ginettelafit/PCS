############################################################
##########         R Code by Ginette Lafit        ##########
##########         ginette.lafit@kuleuven.be      ##########
############################################################

############################################################
############################################################
############################################################
##### PCS k-fold cross-validation for Glasso
############################################################
############################################################
############################################################

# k-folds CV for PCS: Glasso

# Function that estimates the threshold parameter in Partial Correlation Screening when the GMM is estimated using Glasso
# the threshold is selected using k-fold cross-validation
# the input is the (n x p) data matrix (i.e., x)
# n is the number of rows of x
# p is the number of columns of x
# lambda.opt is the regularization parameter to estimate the GGM using Glasso
# fold is the number of folds

# The function returns the threhsold parameter that minimizes the mean squared error

cv_PCS_glasso = function(x,lambda.opt,fold){

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
      sel2 =nn[ !(nn %in% sel) ]
      trainMat[,j] = ind[sel2]
    }
    return(list(trainMat=trainMat,testMat=testMat))
  }

 loss.cv = function(x.train,x.test,lambda.opt,taulist){
  Omegalist = as.matrix(huge(x.train,lambda = lambda.opt ,scr = F, method = "glasso")$icov[[1]])
  R = -cov2cor(Omegalist)
  diag(R) = 1 
  eps = lapply(1:length(taulist),function(i) Psi_Screen_Beta(x.train,R,taulist[i])[[1]])
  loss.re = unlist(lapply(1:length(taulist), function(i) sum(colSums((x.test - x.test%*%eps[[i]])^2))))
  return(loss.re)
 }
  
  n = nrow(x)
  p = ncol(x)
  part.list = cv.part(n, fold)

  taulist = seq(0.0001,1,length=100)

  loss.list = lapply(1:fold, function(k) 
  loss.cv(x[part.list$trainMat[, k], ],x[part.list$testMat[, k], ],lambda.opt,taulist))

  loss.re = matrix(unlist(loss.list), ncol = fold, byrow = FALSE)

  loss.mean = apply(loss.re, 1, mean)
  ind = which.min(loss.mean)
  tau_opt = taulist[ind]

  res = tau_opt
  return(res)  
}