####################################
### Generation of simulated data ###
####################################

DG = function(X,Z,gamma_X,gamma_Z,Sigma_e,outcome="continuous") {

  #library(MASS)
  n = dim(X)[1]
  px = length(gamma_X)
  pz = length(gamma_Z)

  ## treatment
  g_link = X%*%gamma_X + Z%*%gamma_Z
  e = 1/(1+exp(-g_link))         # pass through an inv-logit function
  #e =  1 - exp(-exp(g_link))     # complement log-log model form
  T = rbinom(n,1,e)

  ## outcome
  #### M1
  if(outcome=="continuous") {
    Y = T + X%*%gamma_X + Z%*%gamma_Z + rnorm(n,0,1) }
  #######

  #### M2
  if(outcome=="binary") {
    #T = rbinom(n,1,e)
    pi = 1/(1+exp(-(T+g_link)))
    Y = rbinom(n,1,pi)
  }
  #######

  ## error-prone
  mu_X = rep(0,px)
  err = mvrnorm(n, mu_X, Sigma_e , tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  W = X + err

  ## dataset
  data = cbind(Y,T,W,Z)

  return(data)
}


