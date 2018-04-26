/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.
#include <testthat.h>
#include "SeqSupport.h"


context("Linear Regression") {


  const size_t n=3;
  size_t p=2;
  size_t g=4;
  Eigen::MatrixXd x(n,p);
  Eigen::MatrixXd y(n,g);
  
  x <<0.637091983931509, -2.04219082510305,
    1.38937121784883, 0.0194671703245328,
    1.54582167051811, -0.742589965852777;
  
  y << 0.170885455387323, -1.02734039113708, -0.138839413691922, -0.293390177962924,
    0.712875034397334, 1.18237499560882, -0.127520763981914, 0.92967758762599,
    0.747188638327109, -0.788955024890618, -1.23802307425181, -0.365048855441338;
  
  center_columns(x);
  center_columns(y);
  test_that("column centering works"){
    
    for(int i=0;i<p;i++){
      expect_true(std::fabs(x.col(i).mean())<0.0001);
    }
    for(int i=0;i<g;i++){
      expect_true(std::fabs(y.col(i).mean())<0.0001);
    }
  }

  Eigen::MatrixXd true_se(p,g);
  true_se << 0.051007832056321, 1.58149134877262, 0.724005907496506, 
    1.0078098774105, 0.0908716564928982, 0.447113120031327, 0.428436423698243, 
    0.326323451060123 ;
  Eigen::MatrixXd true_bh(p,g);
  true_bh << 0.661492743490991, 1.1089860138709, -0.821983117570159, 
    0.460679270810818, 0.282235650961964, 0.97662962284506, -0.0857189776936776, 
    0.523774337431585;
  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  Eigen::ArrayXd sx2(p);
  Eigen::MatrixXd BH;
  Eigen::MatrixXd se;
  test_that("mapping bh and se work"){

    calc_sx2(x,sx2);
    calc_BH(x,y,sx2,BH);
    calc_se(x,y,sx2,BH,se);
    for(int i=0; i<p;i++){
      for(int j=0; j<g;j++){
	expect_true(std::fabs(BH(i,j)-true_bh(i,j))<0.00001);
	expect_true(std::fabs(se(i,j)-true_se(i,j))<0.00001);
      }
    }
  }
  p=4;
  g=2;
  Eigen::MatrixXd x2(n,p);
  Eigen::MatrixXd y2(n,g);

  
  x2<< 0.105579390183666, 0.236708748477791, 0.689793618629302, -0.5153607095103, 
    0.893490610034201, -2.12729159298048, -0.327665088031768, 0.37635593215243, 
    -0.65558341230427, 0.995946611587856, -0.191337254834271, 1.34641798251883;
  
  y2<<0.455595824821646, -0.333442793507655, 1.71327440282577, 0.111929799061936, 
    -1.04392085682532, 1.1153636356724;

  center_columns(x2);
  center_columns(y2);


  Eigen::MatrixXd true_bh2(p,g);
  true_bh2<< 1.77882359849931, -0.6406619034298, -0.799206931278637, 
    0.199822815716473, -0.182629326386394, -0.87036449942693, -0.835177367515656, 
    0.781927009930042;
  
  Eigen::MatrixXd true_se2(p,g);
  true_se2<<0.0762744301785894, 0.503714710475512, 0.199255903670028, 
    0.289509641958267, 1.76255455726007, 0.723838604989212, 0.866028147146526, 
    0.108937812938828;
  
  test_that("mapping bh and se work after re-use"){
    calc_sx2(x2,sx2);
    calc_BH(x2,y2,sx2,BH);
    calc_se(x2,y2,sx2,BH,se);
    for(int i=0; i<p;i++){
      for(int j=0; j<g;j++){
	expect_true(std::fabs(BH(i,j)-true_bh2(i,j))<0.00001);
	expect_true(std::fabs(se(i,j)-true_se2(i,j))<0.00001);
      }
    }
  }
    Rcpp::Rcerr<<"Finished C++ tests"<<std::endl;
}

context("Int Ranges") {

  using namespace ranges;

  std::vector<int> el ={1,2,3,4,5,6,7};
  Rcpp::IntegerVector Rel(el.begin(),el.end());

  test_that("We can convert Int Vectors to ranges"){
    auto r = IntegerVector_range(Rel);
    expect_true(r[0]==1);
    expect_true(r.size()==7);
  }

  test_that("We can transform Int Vector ranges"){
    std::vector<int> r = IntegerVector_range(Rel) | view::transform([](int i){
	return(i-1);
      });
    expect_true(r[0]==0);
    expect_true(r[1]==1);
    expect_true(r.size()==7);
  }
  test_that("We can chunk Int Vector ranges"){
    std::vector<int> r = IntegerVector_range(Rel) | view::chunk(2) | view::transform([](auto i){
	return(i.front());
      });
    expect_true(r[0]==1);
    expect_true(r[1]==3);
  }



    Rcpp::Rcerr<<"Finished Ranges tests"<<std::endl;
}
