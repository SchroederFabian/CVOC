#' Expected Loss of the Bayesian Classifier.
#' 
#' @description The function offers a method to select variables by univariate filtering based on the estimated loss of the 
#' univariate Bayesian Classifer. The statistic requires the parametric assumption that the variable consists of a mixture of 
#' Gaussian variables.
#' 
#' @param class a factor vector indicating the class membership of the instances. Must have exactly two levels.
#' @param data a data frame with the variables to filer in columns.
#' @param oc a vector containing three elements. oc[1], the cost of misclassifying a negative instance, oc[2], the cost 
#' of missclassifying a positive instance, and oc[3], the share of negative instances in the population.
#' @param positive a character object indicating the factor label of the positive class.
#' @param robust a logical indicating whether a robust estimator of the mean and variance of the two classes should be used.
#' @param p.val a logical indicating whether the p-values of ebc values under the null hypothesis that both classes are equal
#'  should be calculated. Currently the null distribution is calculated by permutation.
#' @param adj.method a character string indicating the method with which to correct the p-values for multiple testing. See 
#' ?p.adjust.
#' 
#' @return The output is a list with three items. \code{ebc} is a numerical vector containing the ebc values for every 
#' variable of data. \code{p.val} is a numerical vector containing the ebc values for every variable of data. and \code{p.val.adj} is 
#' the corresponding adjusted p-values of ebc. (optional)
#' 
#' @examples 
#' oc <- c(1,3,0.5)
#' class <- factor(c(rep(0, 25), rep(1, 25)), labels=c("neg", "pos"))
#' data <- data.frame("var1"=c(rnorm(25, 0, 1/2), rnorm(25, 1, 2)))
#' res <- ebc(class, data, oc, positive="pos", p.val=TRUE)

#' @export
ebc <- function(class, data, oc=c(1,1,0.5), positive=levels(class)[1], robust=FALSE,  p.val=FALSE, adj.method="BH") {
  
  p <- ifelse(class(data)=="numeric", 1, ncol(data))
  n <- length(class)
  pos <- class==positive
  neg <- !pos
  
  if (nlevels(factor(class))!=2) stop("class must consist of two groups")
  if (class(data)=="numeric") data <- as.data.frame(data)
  if (dim(data)[1] != length(class)) stop("class and data do not correspond.")
  if (!positive%in%levels(class)) stop("positive class doesn't exist.")
  if (length(oc)!=3 | oc[3]<0 | oc[3]>1 | oc[1]<0 | oc[2]<0) stop("oc does not meet requirements.")
  
  epe <- vector(mode="numeric", length=p)
  names(epe) <- colnames(data)
  p.value <- vector(length=length(epe))
  
  for (i in 1:p) {
    
    mu1 <- ifelse(robust, median(data[pos,i], na.rm=TRUE), mean(data[pos,i], na.rm=TRUE))
    mu0 <- ifelse(robust, median(data[neg,i], na.rm=TRUE), mean(data[neg,i], na.rm=TRUE))
    sigma1 <- ifelse(robust, mad(data[pos,i], na.rm=TRUE), sd(data[pos,i], na.rm=TRUE))
    sigma0 <- ifelse(robust, mad(data[neg,i], na.rm=TRUE), sd(data[neg,i], na.rm=TRUE))
    
    if (sigma1==0) sigma1 <- 0.00000000000001
    if (sigma0==0) sigma0 <- 0.00000000000001
    
    if (sigma0==sigma1){
      sigma <- ifelse(robust, mad(data[,i], na.rm=TRUE), sd(data[,i], na.rm=TRUE))
      x <- -( (mu0^2-mu1^2)/(2*sigma^2) - log(oc[1]/oc[2]*oc[3]/(1-oc[3])) ) / ((mu1-mu0)/sigma^2)  
      epe[i] <- ifelse( mu0 < mu1, 	oc[3]*oc[1]*(1-pnorm(x, mu0, sigma)) + (1-oc[3])*oc[2]*pnorm(x, mu1, sigma), 
                        oc[3]*oc[1]*pnorm(x, mu0, sigma)+ (1-oc[3])*oc[2]*(1-pnorm(x, mu1, sigma)))
    } else {
      if (mu0 < mu1) {
        xm <- roots(mu0, mu1, sigma0, sigma1, oc[1], oc[2], oc[3])[2]
        xa <- roots(mu0, mu1, sigma0, sigma1, oc[1], oc[2], oc[3])[1]
        if (is.na(xm)) { epe[i] <- (min((1-oc[3])*oc[2], oc[3]*oc[1])) 
        } else {
          epe[i] <- ifelse (sigma0 < sigma1, 	oc[3] * oc[1] * (1-pnorm(xm, mu0, sigma0) + pnorm(xa, mu0, sigma0))+ (1-oc[3]) * oc[2] * (pnorm(xm, mu1, sigma1)- pnorm(xa, mu1, sigma1)), 
                            oc[3] * oc[1] * (pnorm(xa, mu0, sigma0) - pnorm(xm, mu0, sigma0)) + (1-oc[3]) * oc[2] * (1-pnorm(xa, mu1, sigma1) + pnorm(xm, mu1, sigma1))) 
        }	
      } else {
        xm <- roots(mu0, mu1, sigma0, sigma1, oc[1], oc[2], oc[3])[1]
        xa <- roots(mu0, mu1, sigma0, sigma1, oc[1], oc[2], oc[3])[2]
        if (is.na(xm)) { epe[i] <- (min((1-oc[3])*oc[2], oc[3]*oc[1])) 
        } else {
          epe[i] <- ifelse (sigma0 < sigma1, 	oc[3] * oc[1] * (1-pnorm(xa, mu0, sigma0) + pnorm(xm, mu0, sigma0))+ (1-oc[3]) * oc[2] * (pnorm(xa, mu1, sigma1)- pnorm(xm, mu1, sigma1)), 
                            oc[3] * oc[1] * (pnorm(xm, mu0, sigma0) - pnorm(xa, mu0, sigma0)) + (1-oc[3]) * oc[2] * (1-pnorm(xm, mu1, sigma1) + pnorm(xa, mu1, sigma1)))
        }
      }
      if (sigma0==0 | sigma1==0) {epe[i] <- NA}
      if (is.na(xm)) { epe[i] <- (min((1-oc[3])*oc[2], oc[3]*oc[1]))}
    }
  }
  
  # generate p-values
  if (p.val) {
    reps <- 100000
    data.h0 <- matrix(ncol=reps, nrow=length(class))
    for (j in 1:reps) {data.h0[,j] <- rnorm(length(class))}
    distr.h0 <- ebc(class, data.h0, oc, positive=positive)$ebc
    for (j in 1:ncol(data)) p.value[j] <- sum(distr.h0 <= epe[j])/reps
  } else { p.value <- NULL }
  
  return(list("ebc"=epe, "p.value"=p.value, "adj.p.value"=p.adjust(p.value, method=adj.method)))
}
