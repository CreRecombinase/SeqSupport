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
// #ifdef USE_MKL
// #include "mkl.h"
// #include "mkl_lapacke.h"
// #endif



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


void crossprod_h5(Rcpp::StringVector filenames,Rcpp::StringVector groupnames,Rcpp::StringVector datanames){
  const std::string infile_1=Rcpp::as<std::string>(filenames[0]);
  const std::string infile_2=Rcpp::as<std::string>(filenames[1]);
  const std::string outfile=Rcpp::as<std::string>(filenames[2]);

  const std::string ingroup_1=Rcpp::as<std::string>(groupnames[0]);
  const std::string ingroup_2=Rcpp::as<std::string>(groupnames[1]);
  const std::string outgroup=Rcpp::as<std::string>(groupnames[2]);

  const std::string indata_1=Rcpp::as<std::string>(datanames[0]);
  const std::string indata_2=Rcpp::as<std::string>(datanames[1]);
  const std::string outdata=Rcpp::as<std::string>(datanames[2]);


  using namespace HighFive;
  File file_1(infile_1,File::ReadOnly);
  File file_2(infile_2,File::ReadOnly);

  File out_file(outfile,File::ReadWrite|File::Create);
  int r = register_blosc(nullptr, nullptr);

  auto data_1 = file_1.getGroup(ingroup_1).getDataSet(indata_1);
  auto data_2 = file_2.getGroup(ingroup_2).getDataSet(indata_2);

  auto dims_1 = data_1.getDataDimensions();
  auto dims_2 = data_2.getDataDimensions();

  const size_t rows_1=dims_1[0];
  const size_t cols_1=dims_1[1];

  const size_t rows_2=dims_2[0];
  const size_t cols_2=dims_2[1];

  if(rows_1!=rows_2){
    Rcpp::Rcerr<<"data_1: ("<<ingroup_1<<"/"<<indata_1<<" is "<<rows_1<<" x "<<cols_1<<std::endl;
    Rcpp::Rcerr<<"data_2: ("<<ingroup_2<<"/"<<indata_2<<" is "<<rows_2<<" x "<<cols_2<<std::endl;

    Rcpp::stop("Rows of first dataset must be equal to rows of second dataset");
  }

  const size_t out_rows=cols_1;
  const size_t out_cols=cols_2;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> retmat(out_rows,out_cols);
  Filter filter({out_rows,out_cols},retmat,FILTER_BLOSC, false);
  auto data_out = out_file.createOrGetGroup(outgroup).createDataSet(outdata, DataSpace::From(retmat), AtomicType<double>(), filter.getId(), false);
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> mat_1;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> mat_2;

  data_1.read(mat_1);
  data_2.read(mat_2);

  retmat=mat_1.transpose()*mat_2;
  data_out.write(retmat);

}





//[[Rcpp::export]]
void crossprod_quh_h5(const Rcpp::DataFrame q_dff ,const Rcpp::DataFrame uh_dff,const Rcpp::DataFrame quh_dff){
  register_blosc(nullptr,nullptr);

  const size_t in1_size = q_dff.rows();
  const size_t in2_size = uh_dff.rows();
  const size_t out_size = quh_dff.rows();

  std::unordered_map<std::string,std::shared_ptr<HighFive::File> >  m_file_map;
  std::unordered_map<std::string,std::shared_ptr<HighFive::Group> >  m_group_map;
  std::unordered_map<std::string,std::shared_ptr<HighFive::DataSet> > m_dataset_map;

  if(in1_size!=in2_size){
    Rcpp::stop("input_1 dataframe and input_2 dataframe must have the same number of rows!");
  }
  if(in2_size!=out_size){
    Rcpp::stop("input_1 dataframe and input_2 dataframe must have the same number of rows as output dataframe!");
  }

  MatSlices Q_f(q_dff,m_file_map,m_group_map,m_dataset_map,true);
  MatSlices uh_f(uh_dff,m_file_map,m_group_map,m_dataset_map,true);
  MatSlices Quh_f(quh_dff,m_file_map,m_group_map,m_dataset_map,false);

  using Rowmat= Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
  Rowmat Q;
  Rowmat uh;
  Rowmat Quh;
  Rcpp::Rcerr<<"Using "<<Eigen::nbThreads( )<<" threads"<<std::endl;
  Progress prog_bar(in1_size, true);
  for(int i=0; i<in1_size; i++){
    //    Rcpp::Rcerr<<"Reading Q"<<std::endl;
    Q_f.read(i,Q);
    //    Rcpp::Rcerr<<"Reading uh"<<std::endl;
    uh_f.read(i,uh);
    //    Rcpp::Rcerr<<"Computing crossprod(Q,uh)"<<std::endl;
    Quh.noalias() = Q.transpose()*uh;
    //    Rcpp::Rcerr<<"Writing quh"<<std::endl;
    Quh_f.write(i,Quh);
    prog_bar.increment();
  }
}




void read_ld_chunk_mat_h5(const std::string filename,const int ld_chunk,Eigen::MatrixXd &retmat){

  using namespace HighFive;
  File file(filename,File::ReadOnly);
  std::vector<int> ld_vec;
  file.getGroup("SNPinfo").getDataSet("region_id").read(ld_vec);
  const size_t p=ld_vec.size();
  auto dosage_dat =file.getDataSet("dosage");
  const auto dosage_dims=dosage_dat.getDataDimensions();
  const bool SNPfirst = dosage_dims[0]==p;
  const size_t N = SNPfirst ? dosage_dims[1] : dosage_dims[0];
  if(!SNPfirst){
    if(dosage_dims[1]!=p){
      Rcpp::stop("dosage dims do not match region_id length");
    }
  }
  auto first_it = std::find(ld_vec.begin(),ld_vec.end(),ld_chunk);
  auto last_it = std::find_if_not(first_it,ld_vec.end(),[&](int el){return el == ld_chunk;});
  auto reg_range = std::make_pair(first_it,last_it);
  //auto reg_range = std::equal_range(ld_vec.begin(),ld_vec.end(),ld_chunk);
  if(reg_range.first==ld_vec.end()){
    Rcpp::stop("ld_chunk "+std::to_string(ld_chunk)+"not found!");
  }
  const size_t num_elem=reg_range.second-reg_range.first;
  const size_t SNP_offset = reg_range.first-ld_vec.begin();
  const size_t rownum= SNPfirst ?num_elem : N;
  const size_t colnum= SNPfirst ? N : num_elem;
  const size_t offset_f = SNPfirst ? SNP_offset : 0;
  const size_t offset_s = SNPfirst ?  0 : SNP_offset;
  const size_t chunksize_f = SNPfirst ? num_elem : N;
  const size_t chunksize_s = SNPfirst ?  N : num_elem;
  dosage_dat.selectEigen({offset_f,offset_s},{chunksize_f,chunksize_s},{}).read(retmat);
}


void sim_U(const int p, Eigen::VectorXd tsigu, Eigen::MatrixXd &eX){
  const int g=tsigu.size();
  if(g!=eX.cols() || p!=eX.rows()){
    eX.resize(p,g);
    //    Rcpp::stop("X.cols()!=tvaru.size()");
  }
  for(int i=0; i<g;i++){
    for(int j=0;j<p;j++){
      eX(j, i) = R::rnorm(0,tsigu(i));
    }
  }
}


//[[Rcpp::export(name="sim_U")]]
Eigen::MatrixXd sim_U_exp(const int n, Eigen::VectorXd tsigu){
  const int g=tsigu.size();
  Eigen::MatrixXd X;
  sim_U(n,tsigu,X);
  return(X);
}



  // return(Eigen::EigenMultivariateNormal<double>(Eigen::ArrayXd::Zero(p),tvaru.asDiagonal(),true).samples(n).transpose());
// }

//[[Rcpp::export]]
Eigen::MatrixXd simulate_y_h5(const Rcpp::DataFrame in_dff ,const Rcpp::DataFrame out_dff,const int p, const int N,const int g,Eigen::ArrayXd &tsigu){
  auto r = register_blosc(nullptr,nullptr);

  using Mat=Eigen::MatrixXd;
  std::unordered_map<std::string,std::shared_ptr<HighFive::File> >  m_file_map;
  std::unordered_map<std::string,std::shared_ptr<HighFive::Group> >  m_group_map;
  std::unordered_map<std::string,std::shared_ptr<HighFive::DataSet> > m_dataset_map;
  MatSlices input_f(in_dff,m_file_map,m_group_map,m_dataset_map,true);
  MatSlices output_f(out_dff,m_file_map,m_group_map,m_dataset_map,false);
  Mat X;
  // Mat y(N,g);
  Mat Beta;
  Mat U;
  Mat Y=Mat::Zero(N,g);
  Eigen::VectorXd S;
  const int num_reg=std::distance(input_f.chunk_map.begin(),input_f.chunk_map.end());
  Progress prog_bar(num_reg, true);


  // values near the mean are the most likely
  // standard deviation affects the dispersion of generated values from the mean

  for(auto m_it=input_f.chunk_map.begin();m_it!=input_f.chunk_map.end();m_it++){
    int chunk_id = m_it->first;
    input_f.read_chunk(chunk_id,"dosage",X);
    if(X.rows()!=N){
      X.transposeInPlace();
    }
    X = X.rowwise()-X.colwise().mean();
    S = (1/(X.array().square().colwise().sum()/(N-1)).sqrt())*(1/std::sqrt(N-1));
    output_f.write_chunk(chunk_id,"S",S);
    const size_t tp=X.cols();
    sim_U(tp,tsigu,U);
    Beta=U.array().colwise()*S.array();
    Y=Y+X*Beta;
    Beta.transposeInPlace();
    output_f.write_chunk(chunk_id,"Beta",Beta);
  }
  return(Y);
}


//[[Rcpp::export]]
void map_eQTL_chunk_h5(const Rcpp::DataFrame snp_dff ,const Rcpp::DataFrame exp_dff,const Rcpp::DataFrame uhat_dff,const Rcpp::DataFrame se_dff){

  using Rowmat = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
  using Rowarray = Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;

  using namespace HighFive;
  register_blosc(nullptr,nullptr);

  const size_t snp_rsize = snp_dff.rows();
  const size_t exp_rsize = exp_dff.rows();

  const size_t out_wsize = uhat_dff.rows();

  Rcpp::Rcerr<<"Using "<<Eigen::nbThreads( )<<" threads"<<std::endl;
  if(se_dff.rows()!=out_wsize){
    Rcpp::Rcerr<<"snp_dff has "<<snp_rsize<<" rows"<<std::endl;
    Rcpp::Rcerr<<"exp_dff has "<<exp_rsize<<" rows"<<std::endl;
    Rcpp::Rcerr<<"uhat_dff has "<<out_wsize<<" rows"<<std::endl;
    Rcpp::Rcerr<<"se_dff_dff has "<<out_wsize<<" rows"<<std::endl;
    Rcpp::stop("both uhat_dff and se_dff must have the number of rows equal to nrow(snp_dff)*nrow(exp_dff)");
  }

  const size_t tot_size =snp_rsize*exp_rsize;

  if(tot_size!=out_wsize){
    Rcpp::Rcerr<<"snp_dff has "<<snp_rsize<<" rows"<<std::endl;
    Rcpp::Rcerr<<"exp_dff has "<<exp_rsize<<" rows"<<std::endl;
    Rcpp::Rcerr<<"uhat_dff has "<<out_wsize<<" rows"<<std::endl;
    Rcpp::Rcerr<<"se_dff_dff has "<<se_dff.rows()<<" rows"<<std::endl;

    Rcpp::stop("both uhat_dff and se_dff must have the number of rows equal to nrow(snp_dff)*nrow(exp_dff)");
  }
  std::unordered_map<std::string,std::shared_ptr<HighFive::File> >  m_file_map;
  std::unordered_map<std::string,std::shared_ptr<HighFive::Group> >  m_group_map;
  std::unordered_map<std::string,std::shared_ptr<HighFive::DataSet> > m_dataset_map;

  MatSlices SNP_f(snp_dff,m_file_map,m_group_map,m_dataset_map,true);
  MatSlices EXP_f(exp_dff,m_file_map,m_group_map,m_dataset_map,true);

  std::vector<int> SNP_dims=SNP_f.dims(0);
  std::vector<int> EXP_dims=EXP_f.dims(0);
  MatSlices uh_f(uhat_dff,m_file_map,m_group_map,m_dataset_map,false);
  MatSlices se_f(se_dff,m_file_map,m_group_map,m_dataset_map,false);


  int rk=0;
  Rowmat EXP_chunk;
  Rowarray sy2;
  Rowmat UH_chunk;
  Rowmat se_chunk;
  Rowmat SNP_chunk;
  Eigen::ArrayXd sx2;
  const bool SNP_first= SNP_f.p_first;

  const bool EXP_first= !EXP_f.p_first;

  Progress prog_bar(tot_size, true);
  for(int i=0;i<exp_rsize;i++){
    EXP_f.read(i,EXP_chunk);
    if(EXP_first){
      EXP_chunk.transposeInPlace();
    }
    const int g=EXP_chunk.cols();
    const int N=EXP_chunk.rows();
    EXP_chunk = EXP_chunk.rowwise()-EXP_chunk.colwise().mean();
    sy2=EXP_chunk.array().square().colwise().sum();
    
    for(int j=0;j<snp_rsize;j++){
      SNP_f.read(j,SNP_chunk);
      if(SNP_first){
        SNP_chunk.transposeInPlace();
      }
      const int snp_chunksize=SNP_chunk.cols();
      SNP_chunk = SNP_chunk.rowwise()-SNP_chunk.colwise().mean();

      sx2=SNP_chunk.array().square().colwise().sum();
      UH_chunk=(SNP_chunk.transpose()*EXP_chunk).array().colwise()/sx2;
      se_chunk.resize(snp_chunksize,g);
      if(UH_chunk.rows()!=snp_chunksize || UH_chunk.cols()!=g){
        Rcpp::Rcerr<<"in EXP chunk: "<<i<<std::endl;
        Rcpp::Rcerr<<"in SNP chunk: "<<j<<std::endl;
        Rcpp::Rcerr<<"UH_chunk is :"<<UH_chunk.rows()<<" x "<<UH_chunk.cols()<<std::endl;
        Rcpp::stop("UH_chunk is the wrong dimensions!");

      }
      for(int l=0; l<snp_chunksize;l++){
        for(int k=0; k<g;k++){
          se_chunk(l,k)=std::sqrt((1/(static_cast<double>(N-1)*sx2(l)))*(EXP_chunk.col(k)-(SNP_chunk.col(l)*UH_chunk(l,k))).array().square().sum());
        }
      }
      UH_chunk=UH_chunk.array()/se_chunk.array();
      
      uh_f.write(rk,UH_chunk);
      se_f.write(rk,se_chunk);
      
      rk++;
      
      prog_bar.increment();
    }
  }
}




//[[Rcpp::export]]
Rcpp::NumericMatrix read_ld_chunk_h5(const std::string filename,const int ld_chunk){
  Eigen::MatrixXd retmat;
  read_ld_chunk_mat_h5(filename,ld_chunk,retmat);
  return(Rcpp::wrap(retmat));
}


