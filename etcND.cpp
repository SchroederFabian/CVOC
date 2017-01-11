// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppMP, BH)]]


//#include <stdio.h>
#include <map>
#include <vector>
#include <iterator>
#include <utility>
#include <Rcpp.h>
#include <gmp.h>
//#include <mpfr.h>
#include <RcppMP.h>
#include <gmpxx.h>

using namespace Rcpp;

//[[Rcpp::plugins(cpp11)]]
void addup_fp_left_sub(std::map <std::vector<int>, mpz_class>& memo_fp_left, mpz_class* sum_fp, int level, int start, int stop, int& tp, int& fn, int& tn, int& fp, double& wght) {

  int arr[3] = {level, start, stop};
  std::vector<int> itm(arr, arr+3);
  
  if (memo_fp_left.find(itm) !=memo_fp_left.end()) {
    *sum_fp = *sum_fp + memo_fp_left.at(itm);
  } else {
    
    if (level == 1) {
      *sum_fp = *sum_fp + start - stop + 1;
      memo_fp_left.insert ( std::make_pair(itm, start - stop + 1) );
    } else {

      int level_new = level - 1;
      int mnpsr, mxpsr;
      double crt;
      crt = round((fp - level_new + 1) * wght *1e4)/1e4;
      mnpsr = (fabs(ceil(crt)-crt)<1e-6) ? round(crt+1) : ceil(crt);
      crt = round((tp - fn - mnpsr - (level_new - tn)*wght)*1e4)/1e4;
      mxpsr = floor(crt);
      if (mxpsr >= 0) {
        int start_new = (fp + tp) - mnpsr - (fp - level_new);  
        int stop_new = std::max(level_new, tp - mnpsr + level_new - mxpsr);
        mpz_class bfr(*sum_fp);
        for (int i = start; i >= stop; i--) {
          start_new = std::min(i-1, start_new);
          addup_fp_left_sub(memo_fp_left, sum_fp, level_new, start_new, stop_new, tp, fn, tn, fp, wght);
        }
        mpz_class aftr(*sum_fp);
        memo_fp_left.insert ( std::make_pair(itm, aftr - bfr) );
      } 
      }
    }
  }



//[[Rcpp::plugins(cpp11)]]
void addup_fp_right_sub(std::map <std::vector<int>, mpz_class>& memo_fp_right, mpz_class* sum_fp, int level, int start, int stop, int& tp, int& fn, int& tn, int& fp, double& wght) {
  
  // memoization
  int arr[3] = {level, start, stop};
  std::vector<int> itm(arr, arr+3);
  
  if (memo_fp_right.find(itm) !=memo_fp_right.end()) {
    *sum_fp = *sum_fp + memo_fp_right.at(itm);
  } else {
  
    if (level == 1){
      *sum_fp = *sum_fp + start - stop + 1;
      memo_fp_right.insert ( std::make_pair(itm, start - stop + 1) );
    } else {
      int level_new = level - 1;
      int mnpsr, mxpsr;
      double crt;
      crt = round((fp - level_new + 1) * wght *1e4)/1e4;
      mnpsr = ceil(crt);
      crt = round((tp - fn - mnpsr - (level_new - tn)*wght)*1e4)/1e4;
      mxpsr = (fabs(floor(crt)-crt)<1e-6) ? round(crt - 1) : floor(crt);  
      if (mxpsr >= 0) {
      
        int start_new = (fp + tp) - mnpsr - (fp - level_new);  
        int stop_new = std::max(level_new, tp - mnpsr + level_new - mxpsr);
        mpz_class bfr(*sum_fp);
        for (int i = start; i >= stop; i--) {
          start_new = std::min(i-1, start_new);
          addup_fp_right_sub(memo_fp_right, sum_fp, level_new, start_new, stop_new, tp, fn, tn, fp, wght);
        }
        mpz_class aftr(*sum_fp);
        memo_fp_right.insert ( std::make_pair(itm, aftr - bfr) );
      }  
      }
    }
}

//[[Rcpp::plugins(cpp11)]]
void addup_fn_left_sub(std::map <std::vector<int>, mpz_class>& memo_fn_left, mpz_class* sum_fn, int level, int start, int stop, int& tp, int& fn, int& tn, int& fp, double& wght) {

  // memoization
  int arr[3] = {level, start, stop};
  std::vector<int> itm(arr, arr+3);
  
  if (memo_fn_left.find(itm) !=memo_fn_left.end()) {
    *sum_fn = *sum_fn + memo_fn_left.at(itm);
  } else {
    if (level == 1) {
      *sum_fn = *sum_fn + start - stop + 1;
      memo_fn_left.insert ( std::make_pair(itm, start - stop + 1) );
    } else {
      int level_new = level - 1;
      int mnpsr, mxpsr;
      double crt;
      crt = round((fn - level_new + 1) * wght *1e4)/1e4;
      mnpsr=ceil(crt);
      crt = round( (tn - fp - mnpsr - (level_new - tp) * wght) * 1e4 ) / 1e4;
      mxpsr = floor(crt);
      if (mxpsr >= 0) {
        int start_new = (fn + tn) - mnpsr - (fn - level_new);  
        int stop_new = std::max(level_new, tn - mnpsr + level_new - mxpsr);
        mpz_class bfr(*sum_fn);
        for (int i = start; i >= stop; i--) {
          start_new = std::min(i-1, start_new);
          addup_fn_left_sub(memo_fn_left, sum_fn, level_new, start_new, stop_new, tp, fn, tn, fp, wght);
        } 
        mpz_class aftr(*sum_fn);
        memo_fn_left.insert ( std::make_pair(itm, aftr - bfr) );
      
      } 
      }
    }
}

//[[Rcpp::plugins(cpp11)]]
void addup_fn_right_sub(std::map <std::vector<int>, mpz_class>& memo_fn_right, mpz_class* sum_fn, int level, int start, int stop, int& tp, int& fn, int& tn, int& fp, double& wght) {
  
  // memoization
  int arr[3] = {level, start, stop};
  std::vector<int> itm(arr, arr+3);
  
  if (memo_fn_right.find(itm) !=memo_fn_right.end()) {
    *sum_fn = *sum_fn + memo_fn_right.at(itm);
  } else {
    if (level == 1) {
      *sum_fn = *sum_fn + start - stop + 1;
      memo_fn_right.insert ( std::make_pair(itm, start - stop + 1) );
    } else {

      int level_new = level - 1;
      int mnpsr, mxpsr;
      double crt;
      crt = round((fn - level_new + 1) * wght *1e4)/1e4;
      mnpsr = (fabs(ceil(crt)-crt)<1e-6) ? round(crt+1) : ceil(crt);
      crt = round( (tn - fp - mnpsr - (level_new - tp) * wght) * 1e4 ) / 1e4;
      mxpsr = (fabs(floor(crt)-crt)<1e-6) ? round(crt - 1) : floor(crt);
      if (mxpsr >= 0) {
        int start_new = (fn + tn) - mnpsr - (fn - level_new);  
        int stop_new = std::max(level_new, tn - mnpsr + level_new - mxpsr);
        mpz_class bfr(*sum_fn);
        for (int i = start; i >= stop; i--) {
          start_new = std::min(i-1, start_new);
          addup_fn_right_sub(memo_fn_right, sum_fn, level_new, start_new, stop_new, tp, fn, tn, fp, wght);
        }
        mpz_class aftr(*sum_fn);
        memo_fn_right.insert ( std::make_pair(itm, aftr - bfr) );
      
      }
      }
   }
}


//[[Rcpp::export]]
SEXP addup_fp_left(int level, int start, int stop, int tp, int fn, int tn, int fp, double c0, double c1, double pi0) {
  
  mpz_class sum_fp(0);
  double wght = c0/c1*(tp + fn)/(tn+fp)*pi0/(1-pi0);
  std::map <std::vector<int>, mpz_class> memo_fp_left;
  addup_fp_left_sub(memo_fp_left, &sum_fp, level, start, stop, tp, fn, tn, fp, wght);
  
  //typedef std::map<std::vector<int>, mpz_class>::const_iterator mapIterator;
  //for (mapIterator it = memo_fp_left.begin(); it !=memo_fp_left.end(); ++it)
  //{ Rcpp::Rcout << it->first[1] << ":" << it->first[2] << "::" << it->second << std::endl; }
  
  memo_fp_left.clear();
  return wrap(sum_fp);
}


// [[Rcpp::export]]
SEXP addup_fp_right(int level, int start, int stop, int tp, int fn, int tn, int fp, double c0, double c1, double pi0) {
  
  mpz_class sum_fp(0);
  double wght = c0/c1*(tp+fn)/(tn+fp)*pi0/(1-pi0);
  std::map <std::vector<int>, mpz_class> memo_fp_right;
  addup_fp_right_sub(memo_fp_right, &sum_fp, level, start, stop, tp, fn, tn, fp, wght);
  memo_fp_right.clear();
  return wrap(sum_fp);
}

// [[Rcpp::export]]
SEXP addup_fn_left(int level, int start, int stop, int tp, int fn, int tn, int fp, double c0, double c1, double pi0) {
  
  mpz_class sum_fn(0);
  double wght = c1/c0*(tn+fp)/(tp+fn)*(1-pi0)/pi0;
  std::map <std::vector<int>, mpz_class> memo_fn_left;
  addup_fn_left_sub(memo_fn_left, &sum_fn, level, start, stop, tp, fn, tn, fp, wght);
  memo_fn_left.clear();
  return wrap(sum_fn);
}


// [[Rcpp::export]]
SEXP addup_fn_right(int level, int start, int stop, int tp, int fn, int tn, int fp, double c0, double c1, double pi0) {

  mpz_class sum_fn(0); 
  double wght = c1/c0*(tn+fp)/(tp+fn)*(1-pi0)/pi0;
  std::map <std::vector<int>, mpz_class> memo_fn_right;
  addup_fn_right_sub(memo_fn_right, &sum_fn, level, start, stop, tp, fn, tn, fp, wght);
  memo_fn_right.clear();
  return wrap(sum_fn);
}


