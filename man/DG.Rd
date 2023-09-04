\name{DG}
\alias{DG}
\title{
Generation of artificial data
}
\description{
This function is used to generate an artificial dataset, which contains potential outcomes, treatments, and error-prone confounders.
}
\usage{
DG(X,Z,gamma_X,gamma_Z,Sigma_e,outcome="continuous")
}
\arguments{
  \item{X}{
an \eqn{n\times p_x} matrix of the error-prone confounders
}
  \item{Z}{
an \eqn{n\times p_z} matrix of the precisely measured confounders
}

  \item{gamma_X}{
a \eqn{p_x}-dimensional vector of parameters corresponding to the error-prone confounders X
}
  \item{gamma_Z}{
a \eqn{p_z}-dimensional vector of parameters corresponding to the precisely measured confounders Z
}

  \item{Sigma_e}{
a \eqn{p_x \times p_x} covariance matrix for the classical measurement error model
}

  \item{outcome}{
the indicator of the nature of the outcome variable; \code{outcome="continuous"} reflects normally distributed outcomes; \code{outcome="binary"} gives binary outcomes
}

}
\details{
This function is used to generate artificial data, including potential outcomes, binary treatments, and error-prone and precisely measured confounders.
}

\value{
\item{data}{
an \eqn{n\times (2+p_x+p_z)} matrix of the artificial data. The first column is the potential outcome, and the second column is the binary treatment; column 3 to column (\code{p_x+2}) records error-prone confounders, and the remaining columns record precisely-measured confounders.
}



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
}

\keyword{funtion}
\keyword{core}
