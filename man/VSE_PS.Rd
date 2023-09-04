\name{VSE_PS}
\alias{VSE_PS}
\title{
Variable selection for confounders
}
\description{
This function implements the penalized quadratic loss function to select the informative confounders.
}
\usage{
VSE_PS(V,y,method="lasso",cv="TRUE",alpha=1)
}
\arguments{
    \item{V}{a user-specified matrix in the quadratic loss function
    }

    \item{y}{a vector determined by \link{SIMEX_EST}
    }

    \item{method}{it specifies a choice of the penalty function with options \code{"lasso"} (Tibshirani 1996), \code{"scad"} (Fan and Li 2001) and \code{"mcp"} (Zhang 2010). The default is set as \code{method="lasso"}.
}

    \item{cv}{the usage for choosing the tuning parameter. \code{cv="TRUE"} suggests the use of the cross-validation method, and \code{cv="FALSE"} allows the use of the BIC. The default is set as \code{cv="TRUE"}.
}

    \item{alpha}{the constant appearing in the Elastic Net penalty (Zou and Hastie 2005). The default value is 1.
}





}
\details{
This function is used to do variable selection for informative confounders by various choices of penalty functions.
}

\value{
a vector of estimators in the treatment model, where components with zero values represent confounders that are unimportant and need to excluded; components with nonzero values identify important confounders that enter the treatment model.
}




\references{
\enumerate{
    \item Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized likelihood and its oracle
properties. \emph{Journal of the American Statistical Association}, 96, 1348-1360.

    \item Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. \emph{Journal of the Royal Statistical Society, Series B}, 58, 267-288.

    \item Yi, G. Y. and Chen, L.-P. (2023). Estimation of the average treatment effect with variable selection and measurement error simultaneously addressed for potential confounders. \emph{Statistical Methods in Medical Research}, 32, 691-711.

\item Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax concave penalty. \emph{The Annals of Statistics}, 38, 894-942.

\item  Zou, H., and Hastie, T. (2005). Regularization and variable selection via the elastic net. \emph{Journal of the Royal Statistical Society, Series B}, 67, 301-320.

  }
}

\author{
Chen, L.-P. and Yi, G. Y.
}

\seealso{
  \code{\link{SIMEX_EST}}

}


\examples{

library(MASS)
n = 800
p_x = 10      # dimension of parameters
p_z = 10
p = p_x + p_z
gamma_X = c(rep(1,2),rep(0,p_x-2))
gamma_Z = c(rep(1,2),rep(0,p_z-2))
gamma = c(gamma_X, gamma_Z)

mu_X = rep(0,p_x)
mu_Z = rep(0,p_z)

Sigma_X = diag(1,p_x,p_x)
Sigma_Z = diag(1,p_z,p_z)
Sigma_e = diag(0.2,p_x)
X = mvrnorm(n, mu_X, Sigma_X, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
Z = mvrnorm(n, mu_Z, Sigma_Z, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
data = DG(X,Z,gamma_X,gamma_Z,Sigma_e,outcome="continuous")


y = as.vector(SIMEX_EST(data,PS="logistic",Psi = seq(0,2,length=10),p_x=length(gamma_X),
              K=5, Sigma_e=diag(0.2,p_x)))
V = diag(1,length(y),length(y))

VSE_PS(V,y,method="lasso",cv="TRUE")
VSE_PS(V,y,method="scad",cv="TRUE")
VSE_PS(V,y,method="mcp",cv="TRUE")

}

\keyword{funtion}
\keyword{core}
