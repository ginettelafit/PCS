############################################################
############################################################
#####         R Code by Ginette Lafit        ###############
#####         ginette.lafit@kuleuven.be      ###############
############################################################
############################################################

ls()

rm(list=ls())

library(Rlab)

library(psych)

library(MASS)

library(Tlasso)

library(huge)

set.seed(123)

############################################################
############################################################
######               Data Simulation                 #######
############################################################
############################################################

p = 20 # Set dimenensionality
n = 100 # Set sample size
 
#########################################################################
##### Generate the covariance (Sigma) matrix & simulate data N(0,Sigma) 
#########################################################################

# Model 1: : 2 neighbor Chain Graph

Omega = diag(p)

for (i in 1:p){
Omega[i,i-1] = 0.4
Omega[i-1,i] = 0.4
}

Sigma = solve(Omega)

R = -cov2cor(Omega)
diag(R) = 1

save(R, file = "pcor_1_p_20.RData")

# Simulate Data

x = mvrnorm(n, rep(0, p), Sigma)
x = scale(x, center = TRUE, scale = F)
save(x, file = "sample_1_p_20_n_100.RData")

############################################################

# Model 2: : 3 neighbor Chain Graph

Omega = diag(p)

for (i in 1:p){
Omega[i,i-1] = 0.4
Omega[i-1,i] = 0.4
}

for (i in 2:p){
Omega[i,i-2] = 0.2
Omega[i-2,i] = 0.2
}

Sigma = solve(Omega)

R = -cov2cor(Omega)
diag(R) = 1

save(R, file = "pcor_2_p_20.RData")

# Simulate Data

x = mvrnorm(n, rep(0, p), Sigma)
x = scale(x, center = TRUE, scale = F)
save(x, file = "sample_2_p_20_n_100.RData")

############################################################

# Model 3: : 2 nearest-neighbor graph

Omega = NeighborOmega(p, sd = 1, knn = 2, norm.type = 1)

Sigma = solve(Omega)

R = -cov2cor(Omega)
diag(R) = 1

save(R, file = "pcor_3_p_20.RData")

# Simulate Data

x = mvrnorm(n, rep(0, p), Sigma)
x = scale(x, center = TRUE, scale = F)
save(x, file = "sample_3_p_20_n_100.RData")

############################################################

# Model 4: Random graph

Omega = huge.generator(10000, p, graph = "random",prob=0.1)$omega # if p = 20
Omega = huge.generator(10000, p, graph = "random",prob=0.01)$omega # if p = 60
Omega = huge.generator(10000, p, graph = "random",prob=0.001)$omega # if p = 200

Omega = cov2cor(Omega)

for (i in 1:p){
Omega[i,] = ifelse(abs(Omega[i,])<=1e-6,0,Omega[i,])
}

Sigma = solve(Omega)

R = -cov2cor(Omega)
diag(R) = 1

save(R, file = "pcor_4_p_20.RData")

# Simulate Data

x = mvrnorm(n, rep(0, p), Sigma)
x = scale(x, center = TRUE, scale = F)
save(x, file = "sample_4_p_20_n_100.RData")


