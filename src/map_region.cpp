#include "EigenH5.h"
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen,BH)]]
// [[Rcpp::depends(RcppProgress)]]
//[[Rcpp::depends(RcppParallel)]]

#include <progress.hpp>
#include <progress_bar.hpp>
#include "boost/multi_array.hpp"
//[[Rcpp::plugins(cpp17)]]
#include "eigenmvn.h"
#include <iterator>
#include<H5Tpublic.h>



#include <range/v3/core.hpp>
#include <range/v3/numeric/adjacent_difference.hpp>
#include <range/v3/view/group_by.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/action/transform.hpp>
#include <range/v3/action/join.hpp>
#include <range/v3/view/for_each.hpp>
#include <range/v3/view/zip.hpp>
#include <range/v3/view/zip_with.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/map.hpp>

#include <range/v3/view/indices.hpp>
#include <range/v3/view/chunk.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/span.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/all.hpp>
#include <range/v3/view/c_str.hpp>

#include <range/v3/algorithm/mismatch.hpp>


//
// struct dim_sel{
// public:
//   int in_start;
//   int in_stop;
//   int out_start;
//   int out_stop;
//   bool sorted;
//   int chunksize;
//
//   dim_sel(const int in_start_,int in_stop_,const int out_start_,const int out_stop_,const int dimsize):out_start(out_start_),
// 												       out_stop(out_stop_){
//     if(in_stop_<0){
//       in_stop=dimsize-in_stop;
//     }else{
//       in_stop=in_stop_;
//     }
//     if(in_start_<0){
//       in_start=dimsize-in_start;
//     }else{
//       in_start=in_start_;
//     }
//     if(in_start>in_stop){
//       std::swap(in_start,in_stop);
//       sorted=false;
//     }else{
//       sorted=true;
//     }
//     chunksize =	(in_stop-in_start+1);
//     if(chunksize!=(out_stop-out_start+1)){
//       Rcpp::Rcerr<<"chunksize: "<<chunksize<<", (out_stop-out_start+1): "<<(out_stop-out_start+1)<<std::endl;
//       Rcpp::stop("chunksize mismatch");
//     }
//   }
//   dim_sel(const int offset,const int dimsize){
//     in_stop=dimsize-1;
//     if(offset<0){
//       in_start=dimsize-offset;
//     }else{
//       in_start=offset;
//     }
//     chunksize=in_stop-in_start+1;
//     out_start=0;
//     out_stop=chunksize-1;
//     sorted=true;
//   }
//
//   dim_sel(const int offset, const int chunksize_,const int dimsize):chunksize(chunksize_){
//   if(offset<0){
//     in_start=dimsize-offset;
//   }
//   out_start=0;
//   out_stop=chunksize-1;
//   if(in_start+chunksize>=dimsize){
//     Rcpp::Rcerr<<"offset: "<<offset<<" chunksize: "<<chunksize<<" dimsize: "<<dimsize<<std::endl;
//     Rcpp::stop("Illegal selection, offset+chunksize > dimsize");
//   }
//   sorted=true;
// }
// };
//
// std::vector<std::tuple<int,int,int> > chunk_cont(const Rcpp::IntegerVector inp,int chunksize=0){
//    const int out_size=inp.size();
//   if(chunksize==0){
//     chunksize = out_size;
//   }
//   using namespace std::placeholders;
//   using namespace ranges;
//   auto ir = make_iterator_range(inp.begin(), inp.end());
//   // auto dr = make_iterator_range(diff_v.begin(), diff_v.end());
//   auto b_chunk = std::bind(view::chunk,_1,chunksize);
//   std::vector<std::tuple<int,int,int> > ar= view::zip_with([](int i,int j){
//     return(std::make_tuple(i-1,j));
//   },ir,view::ints(0)) | view::group_by([&](std::tuple<int,int> i, std::tuple<int,int> j){
//     return((std::get<0>(i)-std::get<0>(j))==(std::get<1>(i)-std::get<1>(j)));
//   }) | view::transform(b_chunk) | view::join | view::transform([](auto el){
//     auto elr = el.front();
//     int csize= distance(el);
//     return(std::make_tuple(std::get<0>(elr),std::get<1>(elr),csize));
//   });
//   return(ar);
//
// }


// Rcpp::DataFrame cont_diff(const Rcpp::IntegerVector inp,int chunksize=0){
//   using namespace ranges;


//   int tkk=0;
//   const int n_groups = ar.size();
//   Rcpp::IntegerVector chunk_i(n_groups);
//   Rcpp::IntegerVector in_beg(n_groups);
//   Rcpp::IntegerVector out_beg(n_groups);
//   Rcpp::IntegerVector csize(n_groups);
//   for(int i=0; i<n_groups;i++){
//     auto te=ar[i];
//     chunk_i[i]=i;
//     in_beg[i]=std::get<0>(te);
//     out_beg[i]=std::get<1>(te);
//     csize[i]  = std::get<2>(te);
//   }


//   using namespace Rcpp;

//   return(DataFrame::create( _["chunk_id"]=chunk_i,
//                             _["in_offset"]=in_beg,
//                             _["out_offset"]=out_beg,
//                             _["chunksize"]=csize));
// }


//
//
// template<typename It>
// std::vector<dim_sel> find_cont(const It itb, const It ite,const int total_size, int chunksize=0){
//   using namespace Rcpp;
//   using namespace ranges;
//
//
//
//   const int n_elem = ite-itb;
//
//   Rcpp::IntegerVector te(itb,ite);
//   std::transform(itb,ite,te.begin(),[](int f){return f-1;});
//
//   auto ntbb=te.begin();
//   auto ntb=te.begin();
//   auto nte=te.end();
//   auto it = ntb;
//   int tot_dist=0;
//
//   std::vector<dim_sel> sub_ranges;
//     sub_ranges.reserve(n_elem/2);
//   if(it==nte && chunksize==0){
//     sub_ranges.push_back(dim_sel(0,total_size-1,0,total_size-1,total_size));
//     return(sub_ranges);
//   }
//
//   if(chunksize==0){
//     chunksize = n_elem;
//   }
//
//   while(it!=nte){
//     int sf=0;
//
//     it = std::adjacent_find(ntb,nte,[&](int i,int j){
//       sf++;
//       return(((j-i)!=1) && (sf>=chunksize));
//     });
//     int iti = it==nte ? *(it-1) : *(it);
//     int ntb_pos = ntb-ntbb;
//     int reg_size = it==nte ? it-ntb : (it-ntb+1);
//     sub_ranges.push_back(dim_sel(*ntb,iti,tot_dist,tot_dist+reg_size-1,total_size));
//     if(it!=nte){
//       it++;
//     }
//     tot_dist=tot_dist+reg_size;
//     ntb=it;
//   }
//   return(sub_ranges);
// }
//
//
// Rcpp::DataFrame cont_reg(Rcpp::IntegerVector input_rows=Rcpp::IntegerVector::create(),
// 			 Rcpp::IntegerVector input_cols=Rcpp::IntegerVector::create(),
// 			 Rcpp::IntegerVector chunksizes=Rcpp::IntegerVector::create(0,0),
// 			 Rcpp::IntegerVector dimsize=Rcpp::IntegerVector::create(0,0),
// 			 int chunk_group=0){
//
//   if(input_rows.size()==0){
//     if(dimsize(0)==0){
//       Rcpp::stop("input_rows cannot be empty if	dimsize(0) is not specified");
//     }
//     input_rows=Rcpp::IntegerVector(dimsize(0));
//     std::iota(input_rows.begin(),input_rows.end(),1);
//   }
//
//   if(input_cols.size()==0){
//     if(dimsize(1)==0){
//       Rcpp::stop("input_cols cannot be empty if	dimsize(1) is not specified");
//     }
//     input_cols=Rcpp::IntegerVector(dimsize(1));
//     std::iota(input_cols.begin(),input_cols.end(),1);
//   }
//
//   auto ret_rows=chunk_cont(input_rows,chunksizes(0));
//   auto ret_cols=chunk_cont(input_cols,chunksizes(1));
//   //auto ret_rows=find_cont(input_rows.begin(),input_rows.end(),dimsize(0),chunksizes(0));
//   //  auto ret_cols=find_cont(input_cols.begin(),input_cols.end(),dimsize(1),chunksizes(1));
//
//   const int nret_rows=ret_rows.size();
//   const int nret_cols=ret_cols.size();
//   const	int chunk_num =nret_rows*nret_cols;
//   using namespace Rcpp;
//   Rcpp::IntegerVector row_offsets(chunk_num);
//   Rcpp::IntegerVector row_chunksizes(chunk_num);
//   Rcpp::IntegerVector col_offsets(chunk_num);
//   Rcpp::IntegerVector col_chunksizes(chunk_num);
//   Rcpp::IntegerVector chunk_groupv(chunk_num);
//   int tgrp=0;
//   for(int i=0;i<nret_rows;i++){
//     auto trow=ret_rows[i];
//
//     for(int j=0;j<nret_cols;j++){
//       auto tcol=ret_cols[j];
//       row_offsets[tgrp]=std::get<0>(trow);
//       col_offsets[tgrp]=std::get<0>(tcol);
//       row_chunksizes[tgrp]=std::get<2>(trow);
//       col_chunksizes[tgrp]=std::get<2>(tcol);
//       chunk_groupv[tgrp]=chunk_group;
//       tgrp++;
//       chunk_group++;
//     }
//   }
//
//   return(DataFrame::create( _["row_offsets"] = row_offsets,
//                             _["row_chunksizes"] = row_chunksizes,
//                             _["col_offsets"] = col_offsets,
//                             _["col_chunksizes"] = col_chunksizes,
// 			    _["chunk_group"]=chunk_groupv));
// }
