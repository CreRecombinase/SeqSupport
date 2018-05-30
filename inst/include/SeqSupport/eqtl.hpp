#pragma once

#include "RcppEigen.h"



void center_columns_d(Eigen::MatrixXd &mat);


void calc_vx_d(const Eigen::MatrixXd &centered_mat,
	     Eigen::ArrayXd &vx);



void calc_BH_d(const Eigen::MatrixXd &centered_snpmat,
	     const Eigen::MatrixXd &centered_expmat,
	       const Eigen::ArrayXd &sx2,
	       Eigen::MatrixXd &BH_mat);



void calc_CP_d(const Eigen::MatrixXd &centered_snpmat,
	     const Eigen::MatrixXd &centered_expmat,
	      Eigen::MatrixXd &CP_mat);


void calc_se_d(const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &centered_snpmat,
	     const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &centered_expmat,
	     const Eigen::Array<double,Eigen::Dynamic,1> &sx2,
	     const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &CP_mat,
	     Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &se_mat);



void calc_se_d(const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &centered_snpmat,
	     const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &centered_expmat,
	     const Eigen::Array<double,Eigen::Dynamic,1> &sx2,
	     const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &CP_mat,
	     Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &se_mat);

void calc_uh_d(Eigen::MatrixXd &UH_mat,
	       const Eigen::MatrixXd &se_mat);



std::pair<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> ,Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > calc_bh_se_d(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &snpmat,
																  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &expmat);
