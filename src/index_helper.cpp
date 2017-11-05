#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <algorithm>
#include <vector>


template<class ForwardIteratorComp,class ForwardIteratorInd> void match_sorted(
    ForwardIteratorComp query_begin,
    ForwardIteratorComp query_end,
    ForwardIteratorComp target_begin,
    ForwardIteratorComp target_end,
    ForwardIteratorInd index_begin){
  const size_t query_size=query_end-query_begin;
  auto t_target_begin=target_begin;

  for(auto i=query_begin;i!=query_end;i++){
    t_target_begin=std::lower_bound(t_target_begin,target_end,*i);
    *index_begin++=t_target_begin-target_begin;
  }
}


//[[Rcpp::export(name="match_sorted")]]
Rcpp::IntegerVector match_sorted_exp(const Rcpp::IntegerVector &query,const Rcpp::IntegerVector &target){
  const size_t query_size=query.size();
  Rcpp::IntegerVector index_vec(query_size);
  match_sorted(query.begin(),query.end(),target.begin(),target.end(),index_vec.begin());
  return(index_vec);
}
