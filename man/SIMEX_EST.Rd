\name{SIMEX_EST}
\alias{SIMEX_EST}
\title{
Simulation and extrapolation (SIMEX) for the treatment model
}
\description{
This function employs the simulation and extrapolation (SIMEX) method to correct for the measurement error effects for confounders in the treatment model.
}
\usage{
SIMEX_EST(data, PS="logistic", Psi=seq(0,1,length=10),p_x=p, K=200,
          extrapolate="quadratic", Sigma_e,
          replicate = "FALSE", RM = 0)
}
\arguments{
    \item{data}{
    an \eqn{n \times (p+2)} matrix of the data. The first column records the observed outcome, the second column displays the values for the binary treatment, and the remaining columns store the observed measurements for the confounders.
}

    \item{PS}{
    the specification of a link function in the treatment model. \code{logistic} refers to the logistic regression function, \code{probit} reflects the probit model, and \code{cloglog} gives the complementary log-log regression model.
}

    \item{Psi}{
    a user-specified sequence of non-negative values taken from an interval. The default is set as \code{Psi=seq(0,1,length=10)}.
    }

    \item{p_x}{the dimension of the error-prone confounders
    }

    \item{K}{a user-specified positive integer, with the default value set as 200
    }

    \item{extrapolate}{
    the extrapolation function in Step 3. \code{quadratic} reflects the quadratic polynomial function, \code{linear} gives the linear polynomial function, \code{RL} is the rational linear function, and \code{cubic} refers to the cubic polynomial function.    }

    \item{Sigma_e}{
    the covariance matrix for the measurement error model
}

    \item{replicate}{
    the indicator for the availability of repeated measurements in the confounders. \code{replicate = "FALSE"} refers no repeated measurements and \code{replicate = "TRUE"} indicates that repeated measurements exist in the dataset. The default is set as \code{replicate = "FALSE"}.
}

    \item{RM}{a \eqn{p_x}-dimensional user-specified vector with entries being the number of repetitions for each confounder. For example, \code{RM = c(2,2,3)} indicates that three confounders in X have repeated measurements, where the first and second confounders have two repetitions and the third one has three repetitions. The default of \code{RM} is set as the \eqn{p_x}-dimensional zero vector, i.e., \code{RM = rep(0,p_x)}.
}
}
\details{
This function is used to implement the simulation and extrapolation (SIMEX) method to estimate parameters in the treatment model.
}

\value{
a vector of estimators in the treatment model
}




\references{
Yi, G. Y. and Chen, L.-P. (2023). Estimation of the average treatment effect with variable selection and measurement error simultaneously addressed for potential confounders. \emph{Statistical Methods in Medical Research}, 32, 691-711.
}

\author{
Chen, L.-P. and Yi, G. Y.
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


y = SIMEX_EST(data,PS="logistic",Psi = seq(0,2,length=10),p_x =length(gamma_X),
              K=5, Sigma_e=diag(0.2,p_x))
}

\keyword{funtion}
\keyword{core}
