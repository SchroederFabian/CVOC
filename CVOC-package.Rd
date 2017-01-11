\name{CVOC-package}
\alias{CVOC-package}
\alias{CVOC}
\docType{package}
\title{
  The package provides two filter statistics for classification.
}
\description{
  The functions ebc and etc provide two methods for univariate classification that allow for the consideration of operating conditions. These can be used as a filter, a variable selction method. 
}
\author{
Fabian Schroeder
Maintainer: Fabian Schroeder <fabian.schroeder.fl@ait.ac.at>
}
\references{
  This optional section can contain literature or other references for
  background information.
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}
\examples{
  \dontrun{
  # set parameter of simulation

  set.seed(123)
  oc <- c(1,3,0.5)
  mu <- seq(0.05,2,0.1)
  sigma <- 2^seq(-3,3,0.15)
  n0 <- 70
  n1 <- 45
  p <- 1000	# Anzahl der Noise Features

  # generate data
  
  # response
  class <- factor(c(rep(0, n0), rep(1, n1)), labels=c("neg", "pos"))
  # Signal Features
  data <- matrix(ncol=length(mu)*length(sigma)+p, nrow=n0+n1)
  for (i in 1:length(mu)) {
    for (j in 1:length(sigma)) {
      data[,(i-1)*length(sigma) + j] <- c(rnorm(n0, 0, 1/sigma[j]), rnorm(n1, mu[i], sigma[j]))
    }
  }
  # Noise Features
  sf <- length(mu)*length(sigma)
  for (j in 1:p) {
    data[,sf+j] <- rnorm(n0+n1, 0, 2)
  }

  # apply etc and ebc
  
  result.etc <- etc(class, data, oc, positive="pos", plot=FALSE, p.val=TRUE, adj.method="BH")
  result.ebc <- ebc(class, data, oc, positive="pos", robust=FALSE, p.val=FALSE)


 
  }
}
