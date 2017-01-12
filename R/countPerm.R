
countPerm <- function(left, tp, fn, tn, fp, c0, c1, pi0) {
  
  rnd <- 5
  n1 <- tp + fn
  n0 <- tn + fp
  n <- n0 + n1
  sum.fp <- gmp::as.bigz(0)
  sum.fn <- gmp::as.bigz(0)
  
  # number of permutations of false positives
  
  # starting values
  level.init <- fp
  mnpsr <- ifelse(left, min( which( round((c(1:n)*c1/n1*(1-pi0)),rnd) > round(c0/n0*pi0,rnd))),
                  min( which( round((c(1:n)*c1/n1*(1-pi0)),rnd) >= round(0/n0*pi0,rnd))))
  mxpsr <- ifelse(left, suppressWarnings( max( which( round( (fp - level.init + tn)*c0/n0*pi0  +  (tp - mnpsr - c(0:n))*c1/n1*(1-pi0),rnd)  >= round(fp*c0/n0*pi0 + fn*c1/n1*(1-pi0),rnd))-1)),
                  suppressWarnings( max( which( round( (fp - level.init + tn)*c0/n0*pi0  +  (tp - mnpsr - c(0:n))*c1/n1*(1-pi0),rnd)  >  round(fp*c0/n0*pi0 + fn*c1/n1*(1-pi0),rnd))-1)))
  if (!is.finite(mxpsr)) {return(gmp::as.bigz(0))}
  start.init <- (fp + tp) - mnpsr
  stop.init <- max(level.init, tp - mnpsr + level.init - mxpsr)
  
  # initialize
  if (level.init==0) {
    if (tp==0 & left==FALSE) {
      sum.fp <- gmp::as.bigz(0)
    } else {
      sum.fp <- gmp::as.bigz(1)
    }
  } else {
    if (left) { sum.fp <- addup_fp_left(level.init, start.init, stop.init, tp, fn, tn, fp, c0, c1, pi0)
    } else { sum.fp <- addup_fp_right(level.init, start.init, stop.init, tp, fn, tn, fp, c0, c1, pi0)
    }
  }
  
  # number of permutations for false negatives
  
  # staring values
  level.init <- fn
  mnpsr <- ifelse(left, min( which( round( c(1:n)*c0/n0*pi0,rnd) >= round( c1/n1*(1-pi0),rnd))),
                  min( which( round( c(1:n)*c0/n0*pi0,rnd) >  round( c1/n1*(1-pi0),rnd))))
  mxpsr <- ifelse(left, suppressWarnings( max( which( round( (fn - level.init + tp)*c1/n1*(1-pi0)  +  (tn - mnpsr - c(0:n))*c0/n0*pi0,rnd)  >= round(fp*c0/n0*pi0 + fn*c1/n1*(1-pi0),rnd))-1)),
                  suppressWarnings( max( which( round( (fn - level.init + tp)*c1/n1*(1-pi0)  +  (tn - mnpsr - c(0:n))*c0/n0*pi0,rnd)  >  round(fp*c0/n0*pi0 + fn*c1/n1*(1-pi0),rnd))-1)))
  if (!is.finite(mxpsr)) {return(gmp::as.bigz(0))}
  start.init <- (fn + tn) - mnpsr
  stop.init <- max(level.init, tn - mnpsr + level.init - mxpsr)
  
  # initialize
  if (level.init==0) {
    if (tn==0 & left==FALSE) {
      sum.fn <- gmp::as.bigz(0)
    } else {
      sum.fn <- gmp::as.bigz(1)
    }
  } else {
    if (left) { sum.fn <- addup_fn_left(level.init, start.init, stop.init, tp, fn, tn, fp, c0, c1, pi0)
    } else { sum.fn <- addup_fn_right(level.init, start.init, stop.init, tp, fn, tn, fp, c0, c1, pi0)
    }
  }
  
  outpt <- gmp::mul.bigz(sum.fn, sum.fp)
  return(outpt)
}

