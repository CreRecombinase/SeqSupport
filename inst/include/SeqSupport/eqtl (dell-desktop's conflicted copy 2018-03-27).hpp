#pragma once

#include "RcppEigen.h"


template<typename DerivedA>
void center_columns(Eigen::DenseBase<DerivedA> &mat){
  mat=mat.rowwise()-mat.colwise().mean();
}

template<typename DerivedA, typename DerivedB>
void calc_vx(const Eigen::MatrixBase<DerivedA> &centered_mat,
	      Eigen::ArrayBase<DerivedB> &vx){
  const size_t p=centered_snpmat.cols();
  vx=centered_mat.array().square().colwise().sum();
  if(p!=vx.size()){
    Rcpp::Rcerr<<"matrix has "<<p<<" cols, while vx has "<<vx.size()<<" elements!";
    Rcpp::stop("matrix must have the same number of columns as vx has elements");
  }
}



template<typename DerivedA, typename DerivedB,typename DerivedC,typename DerivedD>
void calc_BH(const Eigen::MatrixBase<DerivedA> &centered_snpmat,
	     const Eigen::MatrixBase<DerivedB> &centered_expmat,
	     const Eigen::ArrayBase<DerivedC> &sx2,
	     Eigen::MatrixBase<DerivedD> &BH_mat){
  const size_t N=centered_snpmat.rows();
  const size_t p=centered_snpmat.cols();
  const size_t g=centered_expmat.cols();
  if(centered_expmat.rows()!=N){
    Rcpp::Rcerr<<"SNP matrix has "<<N<<" rows, while trait matrix has "<<centered_expmat.rows()<<std::endl;
    Rcpp::stop("SNP and EXP matrices must have the same number of rows");
  }
  if(p!=sx2.size()){
    Rcpp::Rcerr<<"SNP matrix has "<<p<<" cols, while sx2 has "<<sx2.size()<<" elements!";
    Rcpp::stop("SNP matrix must have the same number of columns as sx2 has elements");
  }
  BH_mat=(centered_snpmat.transpose()*centered_expmat).array().colwise()/sx2;
  if(BH_mat.rows()!=p || BH_mat.cols()!=g){
    Rcpp::Rcerr<<"Expecting matrix of dimension: "<<p<<"x"<<g<<std::endl;
    Rcpp::Rcerr<<"BH_mat is :"<<BH_mat.rows()<<" x "<<BH_mat.cols()<<std::endl;
    Rcpp::stop("BH_mat is the wrong dimensions!");
  }
  
}




template<typename DerivedA, typename DerivedB,typename DerivedC,typename DerivedD>
void calc_CP(const Eigen::MatrixBase<DerivedA> &centered_snpmat,
	     const Eigen::MatrixBase<DerivedB> &centered_expmat,
	     Eigen::MatrixBase<DerivedD> &CP_mat){
  const size_t N=centered_snpmat.rows();
  const size_t p=centered_snpmat.cols();
  const size_t g=centered_expmat.cols();
  if(centered_expmat.rows()!=N){
    Rcpp::Rcerr<<"SNP matrix has "<<N<<" rows, while trait matrix has "<<centered_expmat.rows()<<std::endl;
    Rcpp::stop("SNP and EXP matrices must have the same number of rows");
  }

  CP_mat=(centered_snpmat.transpose()*centered_expmat);
  if(CP_mat.rows()!=p || CP_mat.cols()!=g){
    Rcpp::Rcerr<<"Expecting matrix of dimension: "<<p<<"x"<<g<<std::endl;
    Rcpp::Rcerr<<"CP_mat is :"<<CP_mat.rows()<<" x "<<CP_mat.cols()<<std::endl;
    Rcpp::stop("CP_mat is the wrong dimensions!");
  }
  
}

template<typename T,int Oa>
void calc_se(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Oa> &centered_snpmat,
	     const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Oa> &centered_expmat,
	     const Eigen::Array<T,Eigen::Dynamic,1,Oa> &sx2,
	     const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Oa> &CP_mat,
	     Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Oa> &se_mat){
}


template<typename T,int Oa>
void calc_se(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Oa> &centered_snpmat,
	     const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Oa> &centered_expmat,
	     const Eigen::Array<T,Eigen::Dynamic,1,Oa> &sx2,
	     const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Oa> &CP_mat,
	     Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Oa> &se_mat){

  const size_t N=centered_snpmat.rows();
  const size_t p=centered_snpmat.cols();
  const size_t g=centered_expmat.cols();

  se_mat.resize(p,g);
  if(se_mat.rows()!=p){
    Rcpp::stop("se_mat resize of rows failed");
  }
  if(se_mat.cols()!=g){
    Rcpp::stop("se_mat resize of cols failed");
  }
  if(BH_mat.rows()!=se_mat.rows() || BH_mat.cols()!=se_mat.cols()){
    Rcpp::Rcerr<<"Expecting matrix of dimension: "<<p<<"x"<<g<<std::endl;
    Rcpp::Rcerr<<"BH_mat is :"<<BH_mat.rows()<<" x "<<BH_mat.cols()<<std::endl;
    Rcpp::Rcerr<<"se_mat is :"<<se_mat.rows()<<" x "<<se_mat.cols()<<std::endl;
    Rcpp::stop("BH_mat/se_mat is the wrong dimensions!");
  }
  for(int l=0; l<p;l++){
    for(int k=0; k<g;k++){
      se_mat(l,k)=std::sqrt((1/(static_cast<T>(N-1)*sx2(l)))*(centered_expmat.col(k)-(centered_snpmat.col(l)*BH_mat(l,k))).array().square().sum());
    }
  }
}

template<typename DerivedA, typename DerivedB>
void calc_uh(Eigen::MatrixBase<DerivedA> &UH_mat,
	     const Eigen::MatrixBase<DerivedB> &se_mat){
  UH_mat=UH_mat.array()/se_mat.array();
}



template<typename T,int Options>
std::pair<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Options>,Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Options>> calc_bh_se(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Options> &snpmat,
									     Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Options> &expmat){

  using Mat= typename Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Options>;
  using Array= typename Eigen::Array<T,Eigen::Dynamic,1,Options>;
  
  center_colums(snpmat);
  center_columns(expmat);
  Array sx2;
  calc_sx2(snpmat,sx2);
  std::pair<Mat,Mat> retmats;
  calc_BH(snpmat,expmat,sx2,retmats.first);
  calc_se(snpmat,expmat,sx2,retmats.first,retmats.second);
  return(retmats);
}
  
