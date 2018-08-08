#include "SeqSupport.h"
// [[Rcpp::depends(RcppEigen,BH,RcppProgress,RcppParallel,EigenH5)]]
//[[Rcpp::plugins(cpp17)]]
#include <RcppParallel.h>
#include<H5Tpublic.h>

// [[Rcpp::depends(RcppParallel)]]

//[[Rcpp::export]]
void ld_p_h5(const std::string input_filename,
	     const std::string output_filename,
	     const std::vector<size_t> subset_snp,
	     const std::vector<size_t> snp_id,
	     const std::vector<size_t> subset_ind,
	     const std::vector<double> map,
	     const size_t chunksize=10000){











}
