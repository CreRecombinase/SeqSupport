#include "EigenH5.h"
#include <eigenh5/indexers.hpp>
//[[Rcpp::depends(xtensor)]]
// [[Rcpp::depends(RcppEigen,BH)]]
// [[Rcpp::depends(RcppProgress)]]
//[[Rcpp::depends(RcppParallel)]]

#include <progress.hpp>
#include <progress_bar.hpp>
#include "boost/multi_array.hpp"
//[[Rcpp::plugins(cpp17)]]

#include <iterator>


// void xtensor_crossprod_quh(const Rcpp::List file_l,const bool doTranspose=true,Rcpp::NumericVector allele_flip=Rcpp::NumericVector::create()){
