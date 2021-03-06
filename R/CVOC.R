#' CVOC
#' 
#' @description 
#' The package contains two functions, ebc and etc. These are methods for univariate binary classification that allow for the consideration of operating conditions and can be used as a filter, a variable selction method. 
#' 
#' @examples
#' # set parameter
#' oc <- c(1,3,0.5)
#' mu <- seq(0.05,2,0.1)
#' sigma <- 2^seq(-3,3,0.15)
#' n0 <- 70
#' n1 <- 45
#' p <- 1000	
#' 
#' # generate data
#' class <- factor(c(rep(0, n0), rep(1, n1)), labels=c("neg", "pos"))
#' data <- matrix(ncol=length(mu)*length(sigma)+p, nrow=n0+n1)
#' for (i in 1:length(mu)) {
#'     for (j in 1:length(sigma)) {
#'         data[,(i-1)*length(sigma) + j] <- c(rnorm(n0, 0, 1/sigma[j]), 
#'                                             rnorm(n1, mu[i], sigma[j]))
#'     }
#' }
#' 
#' sf <- length(mu)*length(sigma)
#' for (j in 1:p) {
#'     data[,sf+j] <- rnorm(n0+n1, 0, 2)
#' }
#' 
#' # apply etc and ebc
#' res.etc <- etc(class, data, oc, positive="pos", p.val=TRUE)
#' res.ebc <- ebc(class, data, oc, positive="pos", robust=FALSE)
#' 
#' 
#' @docType package
#' @name CVOC-package
#' @useDynLib CVOC
#' @importFrom Rcpp evalCpp
#' @importFrom stats stepfun pnorm
#' @importFrom gmp chooseZ as.bigz add.bigz mul.bigz
#' @importClassesFrom gmp bigq bigz

NULL