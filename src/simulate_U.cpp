#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(cpp17)]]
#include <highfive/H5DataSet.hpp>
#include <highfive/H5Filter.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Attribute.hpp>
#include <highfive/H5Utility.hpp>
#include <highfive/H5DataType.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5PropertyList.hpp>
#include <highfive/H5FileDriver.hpp>
#include <highfive/H5Object.hpp>
#include <highfive/H5Selection.hpp>
#include <blosc_filter.h>
#include<H5Tpublic.h>

//[[Rcpp::export]]
void map_eQTL_h5(const Rcpp::StringVector SNP_path,
                 const Rcpp::StringVector EXP_path,
                 const Rcpp::StringVector out_path,
                 int SNP_chunksize=20000,
                 int EXP_chunksize=-1){

  using namespace HighFive;
  File snp_f(Rcpp::as<std::string>(SNP_path[0]),File::ReadOnly);
  File exp_f(Rcpp::as<std::string>(EXP_path[0]),File::ReadOnly);

  DataSet SNP_d=snp_f.getGroup(Rcpp::as<std::string>(SNP_path[1])).getDataSet(Rcpp::as<std::string>(SNP_path[2]));
  DataSet EXP_d=exp_f.getGroup(Rcpp::as<std::string>(EXP_path[1])).getDataSet(Rcpp::as<std::string>(EXP_path[2]));

  std::vector<size_t> SNP_dims=SNP_d.getDataDimensions();
  std::vector<size_t> EXP_dims=EXP_d.getDataDimensions();

  bool SNP_first;
  bool EXP_first;
  size_t p,N,g;

  if(SNP_dims[0]==EXP_dims[0]){
    p=SNP_dims[1];
    N=SNP_dims[0];
    g=EXP_dims[1];
    SNP_first=false;
    EXP_first=false;
  }else{
    if(SNP_dims[1]==EXP_dims[0]){
      p=SNP_dims[0];
      N=SNP_dims[1];
      g=EXP_dims[1];
      EXP_first=false;
      SNP_first=true;
    }else{
      if(SNP_dims[0]==EXP_dims[1]){
        p=SNP_dims[1];
        N=SNP_dims[0];
        g=EXP_dims[0];
        EXP_first=true;
        SNP_first=false;
      }else{
        if(SNP_dims[1]==EXP_dims[1]){
          p=SNP_dims[0];
          N=SNP_dims[1];
          g=EXP_dims[0];
          EXP_first=true;
          SNP_first=true;
        }else{
          Rcpp::stop("SNP and Phenotype datasets have different number of individuals");
        }
      }
    }
  }
  std::vector<size_t> out_dims={p,g};
  double *fbool=nullptr;
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > fake_mat(fbool,out_dims[0],out_dims[1]);
  File out_f(Rcpp::as<std::string>(out_path[0]),File::ReadWrite|File::Create);
  auto out_g =out_f.createOrGetGroup(Rcpp::as<std::string>(out_path[1]));
  DataSpace ds = DataSpace::From(fake_mat, false);
  auto plist = H5Pcreate(H5P_DATASET_CREATE);

  int r = register_blosc(nullptr, nullptr);
  if(r<0){
    Rcpp::stop("register blosc failed!");
  }
  Rcpp::Rcout<<"register blosc worked!"<<r<<std::endl;
  // Create a new file using the default property lists.
  Filter filter({1000, 1000}, fake_mat, FILTER_BLOSC, false);
  // Create a dataset with double precision floating points
  plist = filter.getId();
  DataSet beta_d = out_g.createDataSet("uh", ds, AtomicType<double>(), plist, false);
  DataSet se_d = out_g.createDataSet("se", ds, AtomicType<double>(), plist, false);

  if(EXP_chunksize<0){
    EXP_chunksize=g;
  }
  if(SNP_chunksize<0){
    SNP_chunksize=p;
  }

  Eigen::MatrixXd SNP_chunk;
  Eigen::MatrixXd EXP_chunk;
  Eigen::MatrixXd UH_chunk;
  Eigen::MatrixXd se_chunk;
  Eigen::ArrayXd sx2;
  Eigen::ArrayXd sy2;


  for(size_t SNP_idx=0;SNP_idx<p;SNP_idx+=SNP_chunksize){
    size_t SNP_chunk_r=std::min(SNP_idx+SNP_chunksize,p-SNP_idx);
    if(SNP_first){
      SNP_d.selectEigen({SNP_idx,0},{SNP_chunk_r,N},{}).read(SNP_chunk);
      SNP_chunk.transposeInPlace();
    }else{
      SNP_d.selectEigen({0,SNP_idx},{N,SNP_chunk_r},{}).read(SNP_chunk);
    }
    SNP_chunk = SNP_chunk.rowwise()-SNP_chunk.colwise().mean();

    sx2=SNP_chunk.array().square().colwise().sum();
    assert(SNP_chunk.rows()==N);
    assert(SNP_chunk.cols()==SNP_chunk_r);

    for(size_t EXP_idx=0;EXP_idx<g;EXP_idx+=EXP_chunksize){
      size_t EXP_chunk_r=std::min(EXP_idx+EXP_chunksize,g-EXP_idx);
      if(EXP_first){
        EXP_d.selectEigen({EXP_idx,0},{EXP_chunk_r,N},{}).read(EXP_chunk);
        EXP_chunk.transposeInPlace();
      }else{
        EXP_d.selectEigen({0,EXP_idx},{N,EXP_chunk_r},{}).read(EXP_chunk);
      }
      EXP_chunk = EXP_chunk.rowwise()-EXP_chunk.colwise().mean();
      sy2=EXP_chunk.array().square().colwise().sum();
      assert(EXP_chunk.rows()==N);
      assert(EXP_chunk.cols()==EXP_chunk_r);
      se_chunk.resize(SNP_chunk_r,EXP_chunk_r);
      for(int i=0; i<EXP_chunk_r;i++){
        se_chunk.col(i)=(sy2(i)/(static_cast<double>(N-1)*sx2)).sqrt();
      }
      UH_chunk=((SNP_chunk.transpose()*EXP_chunk).array().colwise()/sx2).array()/se_chunk.array();
      beta_d.selectEigen({SNP_idx,EXP_idx},{SNP_chunk_r,EXP_chunk_r},{}).write(UH_chunk);
      se_d.selectEigen({SNP_idx,EXP_idx},{SNP_chunk_r,EXP_chunk_r},{}).write(se_chunk);
    }
  }
}

//[[Rcpp::export]]
void test_loop(){

  std::vector<int> test_vec(10,0);
  for(int i=0; i<10; test_vec[++i]=i){
    Rcpp::Rcout<<"i is: "<<i<<std::endl;
    Rcpp::Rcout<<"test_vec: "<<std::endl;
    for(auto j: test_vec){
      Rcpp::Rcout<<j<<" ";
    }
    Rcpp::Rcout<<std::endl;
  }
}


