#################################################################
##########         R Code by Ginette Lafit        ###############
##########         ginette.lafit@kuleuven.be      ###############
#################################################################

###########################################################
###########################################################
###########################################################
##### SPACE k-fold cross-validation
###########################################################
###########################################################
###########################################################

# SPACE + CV

# Function that performs k-fold cross-validation on SPACE to select the regularization parameter 

# the input is the (n x p) data matrix (i.e., x)
# n is the number of rows of x
# p is the number of columns of x
# fold is the number of folds

# The function returns 2 regularization parameters: alpha_opt and alpha_ls1
# alpha_opt minimizes the mean squared error
# alpha_ls1 is the regularization parameter using the one-standard-error-rule

cv_space = function(x, fold){

library(space)

x = scale(x)
 
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

  loss.cv = function(x.train,x.test,lamlist){
    vareps = lapply(1:length(lamlist),
    function(i) space.joint(x.train, lam1=lamlist[[i]], lam2=0, weight=NULL, iter=3)[[1]])

    eps = lapply(1:length(lamlist),function(i) Psi_Screen_Beta(x.train,vareps[[i]],0)[[1]])

    loss.re = unlist(lapply(1:length(lamlist), function(i) 
    sum(colSums((x.test - x.test%*%eps[[i]])^2))))
  return(loss.re)
  }
  
  n = nrow(x)
  p = ncol(x)
  part.list = cv.part(n, fold)

  l_max = sqrt(n)*qnorm(1-0.0001/(2*p^2))
  l_min = sqrt(n)*qnorm(1-0.9/(2*p^2))

  lamlist = seq(l_min,l_max,length=100)

  loss.list = lapply(1:fold, function(k) 
  loss.cv(x[part.list$trainMat[, k], ],x[part.list$testMat[, k], ],lamlist))

  loss.re = matrix(unlist(loss.list), ncol = fold, byrow = FALSE)

  loss.mean = apply(loss.re, 1, mean)
  std.error = apply(loss.re, 1, sd)/sqrt(fold)
  ind = which.min(loss.mean)

  loss.mean.ls = loss.mean[ind:length(loss.mean)]
  alpha.list.ls = lamlist[ind:length(loss.mean)]

  loss.mean.ls.max = loss.mean.ls[loss.mean.ls <= loss.mean[ind] + std.error[ind]]
  alpha.list.ls.max = alpha.list.ls[loss.mean.ls <= loss.mean[ind] + std.error[ind]]
  ind.ls = length(loss.mean.ls.max)

  alpha_opt = lamlist[ind]
  alpha_ls1 = alpha.list.ls.max[ind.ls]
  
  res = list(alpha_opt=alpha_opt,alpha_ls1=alpha_ls1)
  return(res)
}

