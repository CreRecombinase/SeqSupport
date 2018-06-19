#include "SeqSupport.h"
// [[Rcpp::depends(RcppEigen,BH,RcppProgress,RcppParallel,EigenH5)]]
//[[Rcpp::plugins(cpp17)]]
#include <RcppParallel.h>
#include<H5Tpublic.h>

// [[Rcpp::depends(RcppParallel)]]


template<typename DerivedA>
void center_columns(Eigen::DenseBase<DerivedA> &mat){
  mat=mat.rowwise()-mat.colwise().mean();
}



//[[Rcpp::export(name="center_columns")]]
Eigen::MatrixXd center_columns_exp( Eigen::MatrixXd mat){
  center_columns(mat);
  return(mat);
}







template<typename DerivedA, typename DerivedB>
void calc_sx2(const Eigen::MatrixBase<DerivedA> &centered_snpmat,
	      Eigen::ArrayBase<DerivedB> &sx2){
  const size_t p=centered_snpmat.cols();
  sx2=centered_snpmat.array().square().colwise().sum();
  if(p!=sx2.size()){
    Rcpp::Rcerr<<"SNP matrix has "<<p<<" cols, while sx2 has "<<sx2.size()<<" elements!";
    Rcpp::stop("SNP matrix must have the same number of columns as sx2 has elements");
  }
}

template<typename DerivedA, typename DerivedB>
void calc_sy(const Eigen::MatrixBase<DerivedA> &centered_expmat,
	      Eigen::ArrayBase<DerivedB> &sy){
  const size_t g=centered_expmat.cols();
  const size_t n=centered_expmat.rows();
  sy=(centered_expmat.array().square().colwise().sum()/(n-1)).sqrt();
  //  Rcpp::Rcout<<sy<<std::endl;
  if(g!=sy.size()){
    Rcpp::Rcerr<<"EXP matrix has "<<g<<" cols, while sy has "<<sy.size()<<" elements!";
    Rcpp::stop("EXP matrix must have the same number of columns as sy has elements");
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

template<typename DerivedA, typename DerivedB>
void calc_uh(Eigen::MatrixBase<DerivedA> &UH_mat,
             const Eigen::MatrixBase<DerivedB> &se_mat){
  UH_mat=UH_mat.array()/se_mat.array();
}




template<typename T,int Oa>
void calc_se(const Eigen::Array<T,Eigen::Dynamic,1,Oa> &sx2,
		 const Eigen::Array<T,Eigen::Dynamic,1,Oa> &sy,
	     Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Oa> &se_mat){

  //  const size_t N=centered_snpmat.rows();
  const size_t p=sx2.size();
  const size_t g=sy.size();

  //  se_mat.resize(p,g);

  se_mat=(1/sx2.sqrt()).rowwise().replicate(g).rowwise()*sy.transpose();
  if(se_mat.rows()!=p){
    Rcpp::stop("se_mat resize of rows failed");
  }
  if(se_mat.cols()!=g){
    Rcpp::stop("se_mat resize of cols failed");
  }


}


// template<typename T,int Options>
// class trait_mapper{

//   Mat snp_mat;
//   Mat exp_mat;
//   Array sx2;
//   Array sy;

// };



template<typename T,int Options>
std::pair< Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Options>,Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Options > >
calc_bh_se(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Options> &snpmat,
           Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Options> &expmat){

  using Mat= typename Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Options>;
  using Array= typename Eigen::Array<T,Eigen::Dynamic,1,Options>;


  center_columns(snpmat);
  center_columns(expmat);
  Array sx2;
  Array sy;
  calc_sx2(snpmat,sx2);
  calc_sy(expmat,sy);
  std::pair<Mat,Mat> retmats;
  calc_BH(snpmat,expmat,sx2,retmats.first);
  calc_se(sx2,sy,retmats.second);
  calc_uh(retmats.first,retmats.second);
  return(retmats);
}



//[[Rcpp::export]]
Rcpp::ListOf<Rcpp::NumericMatrix> calc_uh_se_exp(const Eigen::MatrixXd input_snpmat, const Eigen::MatrixXd input_expmat,bool parallel=false){
  Eigen::MatrixXd snp_cp =input_snpmat;
  Eigen::MatrixXd exp_cp =input_expmat;
  Eigen::ArrayXd ta;
  using namespace Rcpp;
  auto ret = calc_bh_se(snp_cp,exp_cp);
  return(Rcpp::List::create(_["UH"]=Rcpp::wrap(ret.first),_["se"]=Rcpp::wrap(ret.second)));

}








//[[Rcpp::export]]
Eigen::MatrixXd evd_rnorm_i(const Eigen::Map<Eigen::MatrixXd> Q, const Eigen::Map<Eigen::VectorXd> s,const Eigen::Map<Eigen::MatrixXd> vm){
  const size_t p=s.size();
  const size_t n=vm.cols();
  if(vm.rows()!=p){
    Rcpp::stop("rows of vm must match length of s");
  }
  return(Q*s.matrix().asDiagonal()*vm);
}





//[[Rcpp::export]]
void crossprod_quh_h5(const Rcpp::List file_l,const bool doTranspose=true){
  register_blosc(nullptr,nullptr);

  //const Rcpp::DataFrame q_dff ,const Rcpp::DataFrame uh_dff,const Rcpp::DataFrame quh_dff
  auto r = register_blosc(nullptr,nullptr);

  FileManager<true> rf(Rcpp::StringVector::create());
  FileManager<false> wf(Rcpp::StringVector::create());
  DataQueue<2,double,true> Q_f(file_l["Q"],rf);
  DataQueue<2,double,true> uh_f(file_l["uh"],rf);
  DataQueue<2,double,false> Quh_f(file_l["quh"],wf);
  const int in_q_size = Q_f.getNumSelections();
  const int in_u_size = uh_f.getNumSelections();
  const int in_uh_size = Quh_f.getNumSelections();


  if( in_q_size!=in_u_size || in_q_size != in_uh_size){
    Rcpp::stop("size of all elements of file_l must be equal!");
  }

  using Mat= Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>;
  Mat Q;
  Mat uh;
  Mat Quh;
  Rcpp::Rcerr<<"Using "<<Eigen::nbThreads( )<<" threads"<<std::endl;
  Progress prog_bar(in_q_size, true);
  for(int i=0; i<in_q_size; i++){
    //    Rcpp::Rcerr<<"Reading Q"<<std::endl;
    Q_f.readMat(i,Q,false);
    //    Rcpp::Rcerr<<"Reading uh"<<std::endl;
    uh_f.readMat(i,uh,false);
    //    Rcpp::Rcerr<<"Computing crossprod(Q,uh)"<<std::endl;
    if(doTranspose){
    Quh.noalias() = Q.transpose()*uh;
    }else{
      Quh.noalias() = Q*uh;
    }
    //    Rcpp::Rcerr<<"Writing quh"<<std::endl;
    Quh_f.writeMat(i,Quh);
    prog_bar.increment();
  }
}

void sim_U(const int p, Eigen::VectorXd tsigu, Eigen::MatrixXd &eX){


  const int g=tsigu.size();

  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> d{0,1};
  auto norm = [&] (double) {return d(gen);};
  eX = Eigen::MatrixXd::NullaryExpr(p,g,norm).array().rowwise()*tsigu.array().transpose();
}




//[[Rcpp::export(name="sim_U")]]
Rcpp::NumericMatrix sim_U_exp(const int n, Eigen::VectorXd tsigu,const int chunksize=2){
  const int g=tsigu.size();
  Eigen::MatrixXd X;
  auto t0 = tbb::tick_count::now();
  sim_U(n,tsigu,X);
  auto t1 = tbb::tick_count::now();
  using namespace Rcpp;
  double ret= (t1-t0).seconds();

  return(wrap(X));
}


//[[Rcpp::export]]
Rcpp::NumericVector calc_af_h5(const Rcpp::List file_l){
  auto r = register_blosc(nullptr,nullptr);

  FileManager<true> rf(Rcpp::StringVector::create());
  DataQueue<2,double,true> SNP_f(file_l,rf);

  const int num_reg = SNP_f.getNumSelections();

  auto dimvec_SNP = SNP_f.get_selection_dims();

  const int N=dimvec_SNP.front().back();
  const double Nd = 2*N;
  int p=1;
  for(int i=0; i<num_reg;i++){
    auto dvSNP = dimvec_SNP[i].front();
    p+=dvSNP;
  }
  using Mat=Eigen::MatrixXd;

  Mat X;
  std::vector<double> retvec;
  retvec.reserve(p);
  Eigen::VectorXd Af;

  Progress prog_bar(num_reg, true);
  for(int i=0; i<num_reg;i++){
    //  int chunk_id = m_it->first;
    SNP_f.readMat(i,X ,false);
    //    X = X.rowwise()-X.colwise().mean();
    Af = X.array().rowwise().sum()/Nd;
    retvec.insert(retvec.end(),Af.data(),Af.data()+Af.size());
    //      (1/(X.array().square().colwise().sum()/(N-1)).sqrt())*(1/std::sqrt(N-1));
    //    S_f.writeVector(i,Af);
    prog_bar.increment();
  }
  return(Rcpp::wrap(retvec));
}





//[[Rcpp::export]]
Eigen::MatrixXd simulate_y_h5(const Rcpp::List file_l ,const int p, const int N,const int g,Eigen::ArrayXd &tsigu,Rcpp::NumericVector Af=Rcpp::NumericVector::create()){

  auto r = register_blosc(nullptr,nullptr);

  FileManager<true> rf(Rcpp::StringVector::create());
  FileManager<false> wf(Rcpp::StringVector::create());
  DataQueue<2,double,true> SNP_f(file_l["SNP"],rf);

  DataQueue<1,double,false> S_f(file_l["S"],wf);
  DataQueue<2,double,false> Beta_f(file_l["Beta"],wf);
  const bool use_AF=Af.size()>0;

  if(use_AF){
    if(Af.size()!=p){
      Rcpp::stop("If specifying allele frequency(Af), length(Af) must be equal to p");
    }
  }

  if((SNP_f.getNumSelections() != Beta_f.getNumSelections() )|| (Beta_f.getNumSelections() != S_f.getNumSelections() )){
    Rcpp::stop("SNP,Beta and S must all have the same number of chunks in 'calc_af_h5'");
  }
  const int num_reg = SNP_f.getNumSelections();

  auto dimvec_SNP = SNP_f.get_selection_dims();
  auto dimvec_Beta = Beta_f.get_selection_dims();
  auto dimvec_S = S_f.get_selection_dims();




  for(int i=0; i<num_reg;i++){
    auto dvSNP = dimvec_SNP[i].front();
    auto dvBeta = dimvec_Beta[i].back();
    auto dvS = dimvec_S[i].front();
    if((dvSNP!=dvBeta) || (dvS !=dvSNP)){
      Rcpp::Rcerr<<"In selection i:"<<std::endl;
      Rcpp::Rcerr<<"SNP selection size (rows): "<<dvSNP<<std::endl;
      Rcpp::Rcerr<<"Beta selection size (rows): "<<dvBeta<<std::endl;
      Rcpp::Rcerr<<"S selection size (rows): "<<dvS<<std::endl;
      Rcpp::stop("Selection subsets must be equal!");
    }
  }






  using Mat=Eigen::MatrixXd;
  Mat X;
  // Mat y(N,g);
  Mat Beta;
  Mat U;
  Mat Y=Mat::Zero(N,g);
  Eigen::VectorXd S;
  // const int num_reg=std::distance(input_f.chunk_map.begin(),input_f.chunk_map.end());
  Progress prog_bar(p, true);
  int cp=0;
  for(int i=0; i<num_reg;i++){
  //  int chunk_id = m_it->first;
    SNP_f.readMat(i,X ,true);
    X = X.rowwise()-X.colwise().mean();
    const size_t tp=X.cols();
    if(use_AF){
      S.resize(tp);
      std::copy_n(&Af[cp],tp,S.data());
      S=2*(S.array()*(1-S.array()));
      S=(1/(S.array().sqrt()*std::sqrt(N-1)));
    }else{
      S = (1/(X.array().square().colwise().sum()/(N-1)).sqrt())*(1/std::sqrt(N-1));
    }
    cp+=tp;
    sim_U(tp,tsigu,U);

    prog_bar.increment(tp);
    Beta=U.array().colwise()*S.array();
    Y+=X*Beta;
    Beta.transposeInPlace();
    Beta_f.writeMat(i,Beta);
    S_f.writeVector(i,S);
  }
  if(cp!=p){
    Rcpp::stop("number of total elements not equal to p");
  }

  return(Y);
}






//[[Rcpp::export]]
Eigen::MatrixXd est_spve_h5(const Rcpp::List file_l ,const int N,Eigen::ArrayXd &sigu,const double rel_D_cutoff=0.01){

  auto tr = register_blosc(nullptr,nullptr);

  FileManager<true> rf(Rcpp::StringVector::create());
  FileManager<false> wf(Rcpp::StringVector::create());
  //  DataQueue<2,double,true> se_f(file_l["se"],rf);
  DataQueue<2,double,true> uh_f(file_l["uh"],rf);
  // DataQueue<2,double,true> R_f(file_l["R"],rf);
  DataQueue<1,double,true> D_f(file_l["D"],rf);
  DataQueue<2,double,true> quh_f(file_l["quh"],rf);
  // DataQueue<2,double,true> q_f(file_l["Q"],rf);


  //  DataQueue<2,double,false> beta_mean_f(file_l["Beta"],wf);

  using Mat = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;
  Mat uh,quh,pve_est;


  Eigen::VectorXd d,rel_d;
  const int g=sigu.size();

  const int tot_size = uh_f.getNumSelections();
  Progress prog_bar(tot_size*g, true);
  pve_est.resize(tot_size,g);
  for(int i=0;i<tot_size;i++){
    //    se_f.readMat(i,se,false);
    uh_f.readMat(i,uh,false);
    quh_f.readMat(i,quh,false);
    // R_f.readMat(i,r,false);
    // q_f.readMat(i,q,false);
    D_f.readVector(i,d);



    quh=(quh.array())/((N+uh.array().square()).array().sqrt());

    const int tp=uh.rows();


    for(int j=0;j<g;j++){
      const double varui=1/(sigu(j)*sigu(j));
      pve_est(i,j)=(quh.col(j).array().square()*(d.array()/((d.array()+varui).square())).array()).sum();

    // (((beta_mean.col(j)*beta_mean.col(j).transpose()).array()*r.array())/(uh.col(j)*uh.col(j).transpose()).array().sqrt()).sum();
    prog_bar.increment();
    }

  }
  return(pve_est);

}




// Eigen::MatrixXd true_spve_h5(const Rcpp::List file_l ,const int N){

//   auto tr = register_blosc(nullptr,nullptr);

//   FileManager<true> rf(Rcpp::StringVector::create());
//   FileManager<false> wf(Rcpp::StringVector::create());
//   DataQueue<2,double,true> se_f(file_l["se"],rf);
//   DataQueue<2,double,true> uh_f(file_l["uh"],rf);
//   DataQueue<2,double,true> R_f(file_l["R"],rf);
//   DataQueue<2,double,true> quh_f(file_l["quh"],rf);
//   DataQueue<2,double,true> q_f(file_l["Q"],rf);
//   DataQueue<2,double,true> beta_mean_f(file_l["Beta"],rf);

//   using Mat = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;
//   Mat uh,se,q,r,quh,beta_mean,pve_est;
//   const size_t g = se_f.get_selection_dims().front().back();




//   const int tot_size = se_f.getNumSelections();
//   Progress prog_bar(tot_size*g, true);
//   pve_est.resize(tot_size,g);
//   for(int i=0;i<tot_size;i++){
//     se_f.readMat(i,se,false);
//     const int g=se.cols();
//     uh_f.readMat(i,uh,false);
//     quh_f.readMat(i,quh,false);
//     R_f.readMat(i,r,false);
//     q_f.readMat(i,q,false);
//     beta_mean_f.readMat(i,beta_mean,true);


//     uh=((uh.array()*se.array()).array().square()+N*se.array().square());

//     for(int j=0;j<g;j++){
//       pve_est(i,j)=(((beta_mean.col(j)*beta_mean.col(j).transpose()).array()*r.array())/(uh.col(j)*uh.col(j).transpose()).array().sqrt()).sum();
//       prog_bar.increment();
//     }

//   }
//   return(pve_est);
// }










//[[Rcpp::export]]
void map_eQTL_chunk_h5(const Rcpp::List snp_dff ,const Rcpp::List exp_dff,const Rcpp::List uhat_dff,const Rcpp::List se_dff,Rcpp::NumericVector Af=Rcpp::NumericVector::create()){

  using Mat = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;
  const bool EXP_first=true;
  const bool SNP_first=true;

  using namespace HighFive;
  register_blosc(nullptr,nullptr);
  const bool use_AF= Af.size()>0;



  FileManager<true> rf(Rcpp::StringVector::create());
  FileManager<false> wf(Rcpp::StringVector::create());

  DataQueue<2,double,true> SNP_f(snp_dff,rf);
  DataQueue<2,double,true> EXP_f(exp_dff,rf);

  DataQueue<2,double,false> uh_f(uhat_dff,wf);
  DataQueue<2,double,false> se_f(se_dff,wf);

  Rcpp::Rcerr<<"Using "<<Eigen::nbThreads( )<<" threads"<<std::endl;
  const size_t snp_rsize = SNP_f.getNumSelections();
  const size_t exp_rsize = EXP_f.getNumSelections();

  if((uh_f.getNumSelections() != se_f.getNumSelections() )|| (uh_f.getNumSelections() != (SNP_f.getNumSelections()*EXP_f.getNumSelections()))){
    const size_t out_wsize = uhat_dff.size();

    Rcpp::stop("uh chunks must equal se chunks, which must equal (snp chunks)*(trait chunks)");
  }



  const int num_reg = SNP_f.getNumSelections();
  auto dimvec_SNP = SNP_f.get_selection_dims();
  auto dimvec_EXP = EXP_f.get_selection_dims();
  std::vector<int>af_offset(num_reg,0);
  int ttp=0;
  for(int i=0; i<num_reg;i++){
    auto dvSNP = dimvec_SNP[i];
    af_offset[i]=ttp;
    ttp+=dvSNP.front();;
    // for(int j=0; j<dimvec_EXP.size(); j++){
    //   auto dvEXP = dimvec_EXP[j];
  }



  int rk=0;
  Mat EXP_chunk;
  Mat UH_chunk;
  Mat se_chunk;
  Mat SNP_chunk;
  const int tot_size = SNP_f.getNumSelections()*EXP_f.getNumSelections();

  Progress prog_bar(tot_size, true);
  for(int i=0;i<exp_rsize;i++){
    EXP_f.readMat(i,EXP_chunk,EXP_first);
    const int g=EXP_chunk.cols();
    const int N=EXP_chunk.rows();
    for(int j=0;j<snp_rsize;j++){
      SNP_f.readMat(j,SNP_chunk ,SNP_first);
      if(SNP_chunk.rows()!=EXP_chunk.rows()){
        Rcpp::Rcerr<<"In SNP_chunk: "<<j<<" SNP_chunk is :"<<SNP_chunk.rows()<<" x "<<SNP_chunk.cols()<<std::endl;
        Rcpp::Rcerr<<"In EXP_chunk: "<<i<<" EXP_chunk is :"<<EXP_chunk.rows()<<" x "<<EXP_chunk.cols()<<std::endl;
        Rcpp::stop("Chunksize mismatch!");
      }
      const int	tp=SNP_chunk.cols();

      auto uh_se= calc_bh_se(SNP_chunk,EXP_chunk);
      uh_f.writeMat(rk,uh_se.first);
      se_f.writeMat(rk,uh_se.second);
      rk++;
      prog_bar.increment();
    }
  }
}
