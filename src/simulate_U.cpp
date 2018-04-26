#include "SeqSupport.h"
// [[Rcpp::depends(RcppEigen,BH,RcppProgress,RcppParallel,EigenH5)]]
//[[Rcpp::plugins(cpp17)]]

#include<H5Tpublic.h>




//[[Rcpp::export]]
Eigen::MatrixXd evd_rnorm_i(const Eigen::Map<Eigen::MatrixXd> Q, const Eigen::Map<Eigen::VectorXd> s,const Eigen::Map<Eigen::MatrixXd> vm){
  const size_t p=s.size();
  const size_t n=vm.cols();
  if(vm.rows()!=p){
    Rcpp::stop("rows of vm must match length of s");
  }
  //Rcpp::NumericMatrix retmat(p,n);
  //Eigen::Map<Eigen::MatrixXd> mQ(&Q(0,0),p,p);
  //Eigen::Map<Eigen::VectorXd> ms(&s(0),p);
  //Eigen::Map<Eigen::MatrixXd> mX(&vm(0,0),p,n);
  return(Q*s.matrix().asDiagonal()*vm);
}





//[[Rcpp::export]]
void crossprod_quh_h5(const Rcpp::DataFrame q_dff ,const Rcpp::DataFrame uh_dff,const Rcpp::DataFrame quh_dff,const bool doTranspose=true){
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
    if(doTranspose){
    Quh.noalias() = Q.transpose()*uh;
    }else{
      Quh.noalias() = Q*uh;
    }
    //    Rcpp::Rcerr<<"Writing quh"<<std::endl;
    Quh_f.write(i,Quh);
    prog_bar.increment();
  }
}







void sim_U(const int p, Eigen::VectorXd tsigu, Eigen::MatrixXd &eX){
  const int g=tsigu.size();
  if(g!=eX.cols() || p!=eX.rows()){
    eX.resize(p,g);
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
void map_eQTL_chunk_h5(const Rcpp::List snp_dff ,const Rcpp::List exp_dff,const Rcpp::List uhat_dff,const Rcpp::List se_dff,const bool EXP_first=true,const bool SNP_first=true){

  using Mat = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;


  using namespace HighFive;
  register_blosc(nullptr,nullptr);

  const size_t snp_rsize = snp_dff.size();
  const size_t exp_rsize = exp_dff.size();

  const size_t out_wsize = uhat_dff.size();

  Rcpp::Rcerr<<"Using "<<Eigen::nbThreads( )<<" threads"<<std::endl;
  if(se_dff.size()!=out_wsize){
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
    Rcpp::Rcerr<<"se_dff_dff has "<<se_dff.size()<<" rows"<<std::endl;

    Rcpp::stop("both uhat_dff and se_dff must have the number of rows equal to nrow(snp_dff)*nrow(exp_dff)");
  }
  FileManager rf(Rcpp::StringVector::create(),true);
  FileManager wf(Rcpp::StringVector::create(),false);

  DataQueue<2,double> SNP_f(snp_dff,true,rf);
  DataQueue<2,double> EXP_f(exp_dff,true,rf);

  // std::vector<int> SNP_dims=SNP_f.dims(0);
  // std::vector<int> EXP_dims=EXP_f.dims(0);
  DataQueue<2,double> uh_f(uhat_dff,false,wf);
  DataQueue<2,double> se_f(se_dff,false,wf);


  int rk=0;
  Mat EXP_chunk;
  Mat UH_chunk;
  Mat se_chunk;
  Mat SNP_chunk;


  Progress prog_bar(tot_size, true);
  for(int i=0;i<exp_rsize;i++){
    EXP_f.readMat(i,EXP_chunk,EXP_first);
    const int g=EXP_chunk.cols();
    const int N=EXP_chunk.rows();
    for(int j=0;j<snp_rsize;j++){
      SNP_f.readMat(j,SNP_chunk ,SNP_first);
      auto uh_se= calc_bh_se(SNP_chunk,EXP_chunk);
      uh_f.writeMat(rk,uh_se.first);
      se_f.writeMat(rk,uh_se.second);
      rk++;
      prog_bar.increment();
    }
  }
}
