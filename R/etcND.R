#' Generate the Null Distribution of the Threshold Classifier.
#' 
#' @description The function offers an algorithmic aproach to generating the null distribution of the etc classifier under the null hypothesis that the distributions of the positive and the negative class are identical.
#' 
#' @param n0 an integer indicating the number of negative instances in the sample.
#' @param n1 an integer indicating the number of positive instances in the sample.
#' @param c0 the cost of misclassifying a negative instance.
#' @param c1 the cost of misclassifying a positive instance.
#' @param pi0 a real number between 0 and 1 indicating the percentage of negative instances in the population.
#' 
#' @return a list containing three components:
#' \item{val}{ a vector with the number of possible values than etc can take.}
#' \item{pos.perm}{ the number of possible permutations for the given number of positives and negatives.}
#' \item{fav.perm}{ a vector with the number of favorable permutations for every value in val.}
#' 
#' @examples 
#' etc.genND(25, 27, 1, 3, 0.5)
#' 
#' @export

etc.genND <- function(n0, n1, c0, c1, pi0) {
  
  # generate list with pairs (fp, fn), for which to calculate the favorable permutations
  mat <- round(outer(seq(0,n0)/n0*c0*pi0,seq(0,n1)/n1*c1*(1-pi0), FUN="+"),4)
  val <- sort(as.vector(mat)[!duplicated(as.vector(mat))])
  val <- val[val<=min(c1*(1-pi0), c0*pi0)]
  clst <- list()
  for (v in val) { clst[[as.character(v)]] <- t(which(mat==v, arr.ind=TRUE)-1) }
  
  # calculate the number of possible permutations
  pos.perm <- gmp::chooseZ(n0+n1, n0)
  res.1 <- gmp::as.bigz(0)
  res.2 <- gmp::as.bigz(0)
  
  # for every pair in clst calculate the number of favorable permutations
  fav.perm <- rep(gmp::as.bigz(0), length(val))
  
  for (i in 1:length(clst)){
    erg <- gmp::as.bigz(0)
    for (j in 1:ncol(clst[[i]])) {
      
      res.1 <- countPerm(TRUE, as.numeric(n1 - clst[[i]][2,j]), as.numeric(clst[[i]][2,j]), as.numeric(n0 - clst[[i]][1,j]), as.numeric(clst[[i]][1,j]), c0, c1, pi0)
      res.2 <- countPerm(FALSE, as.numeric(n1 - clst[[i]][2,j]), as.numeric(clst[[i]][2,j]), as.numeric(n0 - clst[[i]][1,j]), as.numeric(clst[[i]][1,j]), c0, c1, pi0)
      erg <- gmp::add.bigz(erg, gmp::add.bigz(res.1, res.2))
      
    }
    fav.perm[i] <- erg
  }
  return(list("val"=val, "pos.perm"=pos.perm, "fav.perm"=fav.perm))
}  
