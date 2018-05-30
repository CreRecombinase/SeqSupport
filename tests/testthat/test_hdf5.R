context("hdf5 (conversion from gds)")

test_that("We can convert from SeqArray GDS files",{
  gdsf <- system.file("test_19.gds",package="SeqSupport")
  outf <- tempfile()

  SeqSupport::gds2hdf5(gdsf,outf)
  X <- EigenH5::read_matrix_h5(outf,"/","dosage")
  tgds <-SeqArray::seqOpen(gdsf,readonly = T)
  gdsX <-t(seqGetData(tgds,"$dosage"))
  expect_equal(gdsX,X,check.attributes=F)
  seqClose(tgds)
})


test_that("We can simulate traits with specified PVE",{
  gdsf <- system.file("test_19.gds",package="SeqSupport")
  outf <- tempfile()
  SeqSupport::gds2hdf5(gdsf,outf)
  snp_df <- EigenH5::read_df_h5(outf,"SNPinfo")
  dims <- EigenH5::dim_h5(outf,"dosage")
  N <- dims[2]
  pve <- c(0.5,0.9)
  bias <- c(0)
  nreps <- 1
  beta_h5file <- tempfile()
  X <- EigenH5::read_matrix_h5(outf,"/","dosage")
  #tgds <-SeqArray::seqOpen(gdsf,readonly = T)
  good_snp <- apply(X,1,function(x){all(!is.na(x))})
  sub_snp_df <- dplyr::filter(snp_df,good_snp,!is.na(MAF)) %>% dplyr::filter(MAF>0)

  p <- nrow(sub_snp_df)

  tparam_df <- gen_tparamdf_norm(pve,bias,nreps,n = N,p = p) %>% dplyr::mutate(n=N,p=p)
  g <- nrow(tparam_df)
  ymat <- gen_sim_phenotype_h5(sub_snp_df,outf,beta_h5file,tparam_df,AF=numeric())
  expect_equal(colMeans(ymat),rep(0,g))
  betamat <- EigenH5::read_matrix_h5(beta_h5file,"/","Beta")
  yh <- t(X[sub_snp_df$snp_id,])%*%t(betamat)
  expect_equal(apply(yh,2,var)/apply(ymat,2,var),c(0.5,0.9),tolerance=0.1)
})








