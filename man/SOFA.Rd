\name{SOFA}
\alias{SOFA}
\title{Situational Optimization Function Analysis}
\description{
Estimates dispositional trait using the SOFA approach from Franco & Laros (in preparation).
}
\usage{
SOFA(x, y, method="bayes")
}
\arguments{
   \item{x}{A vector of a situational predictor variable.}
   \item{y}{A vector of a behavioral criterion variable.}
   \item{method}{A string with the method to be used for estimating the dispositional true scores. Use "bayes" for the Bayesian model, "gam" for the GAM model, "kernel" for the Kernel model, and "loess" for the Loess model. Defaults to the Bayesian model.}
}
\value{
A list containing the following components:
  \item{disposition}{The estimates of the dispositional trait.}
  \item{full}{The full set of results returned by each method.}
  \item{pred}{The estimates of the behavioral criterion variable.}
If any of the Maximum-Likelihood methods is used (gam, kernell, or loess), then BIC is returned:
  \item{bic}{The BIC fit index for the model.}
If the Bayesian model is used, then DIC is returned:
  \item{dic}{The DIC fit index for the model.}
}
\examples{
### Data generation
## Uncorrelated x, v, and u
mu <- rep(0,3)
set.seed(seed)
rho <- 0
Sigma <- matrix(c(1,0,rho,
                  0,1,0,
                  rho,0,1),
                nrow=3, ncol=3)
z  <- MASS::mvrnorm(n=n, mu=mu, Sigma=Sigma, empirical=T)
Z <- pnorm(z)

x <- qnorm(Z[,1])                                       # Observations
v <- qnorm(Z[,2])                                       # Stochastic error
u <- truncnorm::qtruncnorm(Z[,3],a=0,b=Inf,0,           # Efficiency
                           max(x) - min(x))

## Cost frontier
fy <- 1.5 ^ (x)
y <- fy + v + u

### Estimation
M1 <- SOFA(x,y,method="bayes")
M2 <- SOFA(x,y,method="gam")
M3 <- SOFA(x,y,method="kernel")
M4 <- SOFA(x,y,method="loess")

## Compare results
cor(data.frame("Bayesian"= M1$disposition),
               "GAM"     = M2$disposition),
               "Kernel"  = M3$disposition),
               "Loess"   = M4$disposition))
}
\references{
  Franco, V.R., & Laros, J.A. (in preparation). "An operationalization
  of Lewin's Equation: The situational optimization function analysis".
  Unsubmitted.
}
