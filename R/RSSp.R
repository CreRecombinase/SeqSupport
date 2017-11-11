

gen_quh_chunk_mat <- function(Zmat,evdf,gw_snpi,percent_variance=1){
  ld_snpi <- RcppEigenH5::read_ivec(evdf,"LDinfo","snp_id")
  p <-length(gw_snpi)
  subset_cols <- RcppEigenH5::match_sorted(gw_snpi,target = ld_snpi)+1
  if(length(subset_cols)==length(ld_snpi)){
    Dvec <- RcppEigenH5::read_dvec(evdf,"EVD","D")
    Q <- RcppEigenH5::read_2d_mat_h5(h5file = evdf,groupname = "EVD",dataname = "Q")
  }else{
    tR <- RcppEigenH5::read_2d_index_h5(
      h5file = evdf,
      groupname = "LD",
      dataname = "R",
      indvec = subset_cols)[subset_cols,]
    evdR <- eigen(tR)
    Dvec <- evdR$values
    Q <- evdR$vectors
  }
  csD <-cumsum(Dvec)/p
  stopifnot(all.equal(csD[p],1))
  pm <- min(p,min(which(csD>=(percent_variance-2*.Machine$double.eps))))
  quh <- crossprod(Q[,1:pm],Zmat)
  colnames(quh) <- colnames(Zmat)
  return(list(quh=quh,D=Dvec[1:pm],snp_id=gw_snpi))
}

gen_quh_chunk <- function(gwas_df,evdf){

  ld_snpi <- RcppEigenH5::read_ivec(evdf,"LDinfo","snp_id")
  gw_snpi <- gwas_df$snp_id
  p <-length(gw_snpi)
  subset_cols <- RcppEigenH5::match_sorted(gw_snpi,target = ld_snpi)+1
  if(length(subset_cols)==length(ld_snpi)){
    Dvec <- RcppEigenH5::read_dvec(evdf,"EVD","D")
    Q <- RcppEigenH5::read_2d_mat_h5(h5file = evdf,groupname = "EVD",dataname = "Q")
  }else{
    tR <- RcppEigenH5::read_2d_index_h5(
      h5file = evdf,
      groupname = "LD",
      dataname = "R",
      indvec = subset_cols)[subset_cols,]
    evdR <- eigen(tR)
    Dvec <- evdR$values
    Q <- evdR$vectors
  }
  stopifnot(all.equal(sum(Dvec),p))
  LD_info <- mutate(gwas_df,quh=crossprod(Q,Z),D=Dvec)
  return(LD_info)
}
