#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

//[[Rcpp::export]]
Eigen::MatrixXd haplo_2_geno(const Eigen::MatrixXd haplo, bool snps_in_rows=true){
  const bool sir=snps_in_rows;
  if(sir){
    const size_t n=haplo.cols()/2;
    const size_t p=haplo.rows();
    Eigen::MatrixXd retmat(p,n);
    for(size_t i=0;i<n;i++){
      retmat.col(i)=haplo.block(0,i*2,p,2).rowwise().sum();
    }
    return(retmat);
  }else{
    const size_t n=haplo.rows()/2;
    const size_t p=haplo.cols();
    Eigen::MatrixXd retmat(p,n);
    for(size_t i=0; i<n; i++){
      retmat.col(i)=haplo.block(i*2,0,2,p).colwise().sum();
    }
    return(retmat);
  }
}


//[[Rcpp::export]]
Eigen::MatrixXi haplo_2_geno_i(const Eigen::MatrixXi haplo, bool snps_in_rows=true){
  const bool sir=snps_in_rows;
  if(sir){
    const size_t n=haplo.cols()/2;
    const size_t p=haplo.rows();
    Eigen::MatrixXi retmat(p,n);
    for(size_t i=0;i<n;i++){
      retmat.col(i)=haplo.block(0,i*2,p,2).rowwise().sum();
    }
    return(retmat);
  }else{
    const size_t n=haplo.rows()/2;
    const size_t p=haplo.cols();
    Eigen::MatrixXi retmat(p,n);
    for(size_t i=0; i<n; i++){
      retmat.col(i)=haplo.block(i*2,0,2,p).colwise().sum();
    }
    return(retmat);
  }
}

