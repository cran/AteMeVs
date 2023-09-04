
SIMEX_EST = function(data, PS="logistic", Psi=seq(0,1,length=10),p_x=p, K=200, extrapolate="quadratic", Sigma_e,
                     replicate = "FALSE", RM = 0) {  # Step 3 in SIMEX algorithm

  n = dim(data)[1]
  p = dim(data)[2]-2

  #######  inside function  #######

  SIMEX_S1 = function(data, PS, psi,p_x, Sigma_e, replicate,RM) {   # Step 1 in the SIMEX algorithm

    ## Algorithm settings

    Y = data[,1]
    T = data[,2]

    if(p_x < p && p_x > 0) {
      if(replicate == "FALSE")  {
        W = data[,3:(p_x+2)]
        Z = data[,(p_x+3):(p+2)]
        #set.seed(k)
        mu_X = rep(0,p_x)
        e = mvrnorm(n,  mu_X, Sigma_e, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
        W_S = W + sqrt(psi) * e
        x = cbind(W_S,Z)
      }
      ##################################
      if(replicate == "TRUE")    {
        a = 0
        b = RM[1]
        W_S = NULL
        for(l in 1:length(RM)) {
          reme = c(RM,0)
          W = data[,(3+a):(b+2)]
          dij = matrix(rnorm(n*RM[l],0,1),n,RM[l])
          dibar = rowMeans(dij)
          cij = (dij - dibar) / sqrt( rowSums((dij - dibar)^2))
          work = rowMeans(W) + sqrt(psi/RM[l]) * rowSums(cij * W)
          W_S = cbind(W_S,work)
          a = a + reme[l]
          b = b + reme[l+1]
        }
        x = cbind(W_S,Z)
      }
    }

    if(p_x == p) {
      if(replicate == "FALSE") {
        W = data[,3:(p+2)]
        #set.seed(k)
        mu_X = rep(0,p_x)
        e = mvrnorm(n,  mu_X, Sigma_e, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
        W_S = W + sqrt(psi) * e
        x = W_S
      }
      ##################################
      if(replicate == "TRUE")    {
        a = 0
        b = RM[1]
        W_S = NULL
        for(l in 1:length(RM)) {
          reme = c(RM,0)
          W = data[,(3+a):(b+2)]
          dij = matrix(rnorm(n*RM[l],0,1),n,RM[l])
          dibar = rowMeans(dij)
          cij = (dij - dibar) / sqrt( rowSums((dij - dibar)^2))
          work = rowMeans(W) + sqrt(psi/RM[l]) * rowSums(cij * W)
          W_S = cbind(W_S,work)
          a = a + reme[l]
          b = b + reme[l+1]
        }
        x = cbind(W_S,Z)
      }
    }

    if(p_x == 0) {

      W = data[,3:(p+2)]
      x = W
    }
    ## End of settings


    t = T
    Da = cbind(t,x)
    Da = data.frame(Da)
    if(PS == "logistic") {
      output = glm(t~., data = Da,family = "binomial")$coef[-1]
    }

    if(PS == "probit") {
      output = glm(t~.,data = Da, family=binomial(link="probit"))$coef[-1]
    }

    if(PS == "cloglog") {
      output = glm(t~.,data = Da, family=binomial(link="cloglog"))$coef[-1]
    }

    return(output)
  }

  SIMEX_S2 = function(data,PS,Psi,p_x, K, Sigma_e, replicate,RM) {   # Step 2 in SIMEX algorithm
    GAMMA_Psi = NULL

    for(psi in 1 : length(Psi)) {
      G_coll = NULL                   # collection for all k in the following loop

      for(k in 1 : K) {

        GAMMA_psi_k = SIMEX_S1(data, PS, Psi[psi],p_x, Sigma_e, replicate,RM)

        G_coll = rbind(G_coll, GAMMA_psi_k)

      }
      GAMMA_Psi[[psi]] = G_coll

    }

    return(GAMMA_Psi)

  }

  ######################################

  result = SIMEX_S2(data,PS,Psi,p_x,K, Sigma_e, replicate,RM)

  vt = NULL
  for(psi in 1:length(Psi)) {

    vt = rbind(vt,colMeans(result[[psi]]))

  }

  if(extrapolate=="cubic") {

    B = cbind(1,Psi,Psi^2,Psi^3)

    beta = solve(t(B)%*%B) %*% (t(B)%*%vt)

    ext = c(1,-1,1,-1)

    output = ext %*%beta
  }

  if(extrapolate=="quadratic") {

    B = cbind(1,Psi,Psi^2)

    beta = solve(t(B)%*%B) %*% (t(B)%*%vt)

    ext = c(1,-1,1)

    output = ext %*%beta
  }

  if(extrapolate=="linear") {

    B = cbind(1,Psi)

    beta = solve(t(B)%*%B) %*% (t(B)%*%vt)

    ext = c(1,-1)

    output = ext %*%beta
  }


  if(extrapolate=="RL") {

    output = NULL
    for(i in 1:p) {
      Yresponse = vt[,i]
      RL = function(beta) {

        B1 = beta[1]
        B2 = beta[2]
        B3 = beta[3]

        Q = -mean((Yresponse - B1  - (B2 / (B3 + Psi)) )^2 )

        return (Q)

      }

      betahat = optim(c(1,1,1),RL,gr = NULL, method = "Nelder-Mead",control = list(fnscale = -1))$par

      output = c(output,  (  betahat[1]  + (betahat[2] / (betahat[3] -1))   )     )

    }


  }


  return(output)

}
