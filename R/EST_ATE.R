#################################################
####   SIMEX Algorithm: estimation for ATE   ####
#################################################

EST_ATE = function(data, PS="logistic", Psi=seq(0,1,length=10), K=200, gamma,p_x=p, extrapolate="quadratic",
                   Sigma_e, replicate = "FALSE", RM = 0, bootstrap = 100) {

  n = dim(data)[1]
  p = dim(data)[2]-2

  tau_S1 = function(data, psi, gamma,Sigma_e) {   # Step 1 in the SIMEX algorithm

    ## Algorithm settings

    Y = data[,1]
    T = data[,2]

    if(p_x < p && p_x > 0) {
      if(replicate == "FALSE")            {
        W = data[,3:(p_x+2)]
        Z = data[,(p_x+3):(p+2)]
        #set.seed(k)
        mu_X = rep(0,p_x)
        e = mvrnorm(n,  mu_X, Sigma_e, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
        W_S = W + sqrt(psi) * e
        x = matrix(unlist(cbind(W_S,Z)),n,p)
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
        x = matrix(unlist(cbind(W_S,Z)),n,p)
      }
    }

    if(p_x == p) {
      if(replicate == "FALSE")     {
        W = data[,3:(p+2)]
        #set.seed(k)
        mu_X = rep(0,p_x)
        e = mvrnorm(n,  mu_X, Sigma_e, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
        W_S = W + sqrt(psi) * e
        x = matrix(unlist(W_S),n,p)
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
        x = matrix(unlist(W_S),n,p)
      }
    }

    if(p_x == 0) {

      W = data[,3:(p+2)]
      x = matrix(unlist(W),n,p)
    }

    ## End of settings

    g_link = x%*%gamma

    if(PS == "logistic") {
      ps = 1/(1+exp(-g_link)) }

    if(PS == "probit") {
      ps = pnorm(g_link)
      ps[which(ps==1)] = 0.999  }

    if(PS == "cloglog") {
      ps = 1 - exp(-exp(g_link))
      ps[which(ps==1)] = 0.999  }

    output = (sum(T*Y/ps) / sum(T/ps)) - (sum((1-T)*Y/(1-ps)) / sum((1-T)/(1-ps)) )

    return(output)
  }

  tau_S2 = function(data,Psi,K, gamma,Sigma_e) {   # Step 2 in SIMEX algorithm
    GAMMA_Psi = NULL

    for(psi in 1 : length(Psi)) {
      G_coll = NULL                   # collection for all k in the following loop

      for(k in 1 : length(K)) {

        GAMMA_psi_k = tau_S1(data,Psi[psi], gamma,Sigma_e)
        G_coll = rbind(G_coll, GAMMA_psi_k)

      }
      GAMMA_Psi[[psi]] = G_coll

    }

    return(GAMMA_Psi)

  }

  Output = NULL
  for(boot in 1:bootstrap) {

    result = tau_S2(data,Psi,K, gamma,Sigma_e)

    vt = NULL
    for(psi in 1:length(Psi)) {

      vt = rbind(vt,colMeans(result[[psi]]))

    }

    if(extrapolate=="cubic") {
      B = cbind(1,Psi,Psi^2,Psi^3)
      beta = solve(t(B)%*%B) %*% (t(B)%*%vt)
      ext = c(1,-1,1,-1)
      output = ext %*%beta
      #return(output)
    }

    if(extrapolate=="quadratic") {
      B = cbind(1,Psi,Psi^2)
      beta = solve(t(B)%*%B) %*% (t(B)%*%vt)
      ext = c(1,-1,1)
      output = ext %*%beta
      #return(output)
    }

    if(extrapolate=="linear") {
      B = cbind(1,Psi)
      beta = solve(t(B)%*%B) %*% (t(B)%*%vt)
      ext = c(1,-1)
      output = ext %*%beta
      #return(output)
    }

    if(extrapolate=="RL") {
      RL = function(beta) {

        B1 = beta[1]
        B2 = beta[2]
        B3 = beta[3]

        Q = -mean((vt - B1  - (B2 / (B3 + Psi)) )^2 )

        return (Q)

      }
      betahat = optim(c(2,2,2),RL,gr = NULL, method = "Nelder-Mead",control = list(fnscale = -1))$par
      output = betahat[1]  + (betahat[2] / (betahat[3] -1))
      #return(output)
    }

    Output = c(Output,output)

  }

  if(bootstrap > 1) {

    EST = mean(Output)
    VAR = var(Output)/bootstrap
    pvalue = 2*pnorm(-abs(EST/sqrt(VAR)))

    result = t(as.matrix(c(EST,VAR,pvalue)))
    colnames(result) = c("estimate","variance","p-value")

    return(result)

  }

  if(bootstrap == 1) {

    return(Output)

  }


}



