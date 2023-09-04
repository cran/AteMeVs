
VSE_PS = function(V,y,method="lasso",cv="TRUE",alpha=1) {   # y is SIMEX estimator and V is weighted matrix

  #library(ncvreg)
  #### glmnet package
  #var_sel = glmnet(V, y , family = "gaussian", alpha = 1)
  #lambda = which.max(obj2)       # search max
  #lambda = which.max(obj2[7:length(obj2)])+6   # search max
  #coef = coef[, lambda]
  #obj2 = -dev2 + log(n) * reg.df2   # - 2*logL + 2*log(n) x df  in Step 4 by Yi and Chen (201x)
  #########

  if(cv=="TRUE") {

    if(method =="lasso"){
      opt = cv.ncvreg(V, y , family = "gaussian",penalty="lasso" ,alpha )$lambda
      lambda = which(opt == min(opt))
      var_sel = ncvreg(V, y , family = "gaussian",penalty="lasso" ,alpha )    }

    if(method =="scad"){
      opt = cv.ncvreg(V, y , family = "gaussian",penalty="SCAD" ,alpha, gamma=3.7 )$lambda
      lambda = which(opt == min(opt))
      var_sel = ncvreg(V, y , family = "gaussian",penalty="SCAD" ,alpha, gamma=3.7 )     }

    if(method =="mcp"){
      opt = cv.ncvreg(V, y , family = "gaussian",penalty="MCP" ,alpha, gamma=3 )$lambda
      lambda = which(opt == min(opt))
      var_sel = ncvreg(V, y , family = "gaussian",penalty="MCP" ,alpha, gamma=3 )      }

    coef = rbind(var_sel$a0,as.matrix(var_sel$beta))
    coef = coef[ , lambda]

  }

  if(cv=="FALSE") {
    #### ncvreg package   # for lasso, both results are similar.
    if(method =="lasso"){
      var_sel = ncvreg(V, y , family = "gaussian",penalty="lasso" ,alpha ) }

    if(method =="scad"){
      var_sel = ncvreg(V, y , family = "gaussian",penalty="SCAD" ,alpha, gamma=3.7 )}

    if(method =="mcp"){
      var_sel = ncvreg(V, y , family = "gaussian",penalty="MCP" ,alpha, gamma=3 )}



    coef = rbind(var_sel$a0,as.matrix(var_sel$beta))  # extract coefficients at all values of lambda,  including the intercept
    dev2 = deviance(var_sel)
    reg.df2 = var_sel$df


    obj2 = BIC(var_sel)

    lambda = which.min(obj2)       # search max
    coef = coef[ , lambda]
  }



  output = coef[-1]
  output[which(abs(output)<0.5)]=0
  return(output)


}

