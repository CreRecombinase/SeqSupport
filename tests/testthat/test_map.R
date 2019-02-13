context("mapping")



test_that("QR works (kinda) like R",{
  N <- 1000
  cvrt <-rep(1,N)
  qrr <- qr(cvrt)
  qrq <- qr.Q(qrr)
  attr(cvrt,"dim") <- c(N,1L)
  qrc <- orthogonalize_covars(cvrt)
  expect_equal(qrc,qrq)
})



lm_wrapper <-function(xmat,ymat){
  p <- nrow(xmat)
  g <-nrow(ymat)
  n <- ncol(xmat)
  stopifnot(ncol(ymat)==n)
  st_x <- apply(t(xmat),2,function(x){x-mean(x)})
  st_y <- apply(t(ymat),2,function(x){x-mean(x)})
  bh <- matrix(0,p,g)
  se <- matrix(0,p,g)
  for(i in 1:p){
    for(j in 1:g){
      tl <- coef(summary(lm(st_y[,j]~st_x[,i]+0)))
      bh[i,j] <- tl[1]
      se[i,j] <- tl[2]
    }
  }
  return(list(uh=bh/se,se=se))
}


matrix_eqtl_wrapper <-function(xmat,ymat,cvrt_mat=NULL){

  useModel = MatrixEQTL::modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  tempfile_X <- tempfile()
  tempfile_Y <- tempfile()

  output_file_name = tempfile();
  pvOutputThreshold = 0;
  errorCovariance = numeric();
  N <- ncol(xmat)
  p <- nrow(xmat)
  g <- nrow(ymat)

  dplyr::as_data_frame(xmat) %>%  magrittr::set_colnames(paste0("Sam_",1:N)) %>%
    dplyr::mutate(snpid=1:p) %>%
    dplyr::select(snpid,dplyr::starts_with("Sam")) %>% readr::write_delim(tempfile_X,delim="\t")

  dplyr::as_data_frame(ymat) %>%  magrittr::set_colnames(paste0("Sam_",1:N)) %>%
    dplyr::mutate(geneid=1:g) %>%
    dplyr::select(geneid,dplyr::starts_with("Sam")) %>% readr::write_delim(tempfile_Y,delim="\t")

  if(!is.null(cvrt_mat)){
    ncvrt <- nrow(cvrt_mat)
    tempfile_Z <- tempfile()
    dplyr::as_data_frame(cvrt_mat) %>%  magrittr::set_colnames(paste0("Sam_",1:N)) %>%
      dplyr::mutate(id=1:ncvrt) %>%
      dplyr::select(id,dplyr::starts_with("Sam")) %>% readr::write_delim(tempfile_Z,delim="\t")
  }else{
    tempfile_Z <- character()
  }


  snps = MatrixEQTL::SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(tempfile_X);


  gene = MatrixEQTL::SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(tempfile_Y);

  ## Load covariates

  cvrt = MatrixEQTL::SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(tempfile_Z)>0){
    cvrt$LoadFile(tempfile_Z);
  }



  ## Run the analysis

  me = MatrixEQTL::Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = 1,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);

  me_df <-me$all$eqtls %>%  dplyr::mutate(snps=as.integer(as.character(snps)),gene=as.integer(as.character(gene)))
  uh_mat <- dplyr::arrange(me_df,snps,gene) %>%
    dplyr::select(snps,gene,statistic) %>%
    tidyr::spread(key = gene,value=statistic) %>%
    dplyr::select(-snps) %>%
    data.matrix()
  se_mat <-dplyr::arrange(me_df,snps,gene) %>% dplyr::mutate(se=beta/statistic) %>%
    dplyr::select(snps,gene,se) %>%
    tidyr::spread(key = gene,value=se) %>%
    dplyr::select(-snps) %>%
    data.matrix()
  attr(uh_mat,"dimnames") <- NULL
  attr(se_mat,"dimnames") <- NULL
  return(list(uh=uh_mat,se=se_mat))
}

eqtl_wrapper <-function(xmat,ymat,cvrt_mat=matrix(1,ncol(xmat),1),snp_chunks=2,exp_chunks=2,SNPfirst=TRUE,EXPfirst=TRUE){
  tempf_X <- tempfile()
  tempf_Y <- tempfile()
  tempf_Z <- tempfile()
  uhf <- tempf_Z
  p <- nrow(xmat)
  n <- ncol(xmat)
  g <- nrow(ymat)


  stopifnot(ncol(ymat)==n,nrow(cvrt_mat)==n)

  if(!SNPfirst){
    xmat <- t(xmat)
  }
  if(!EXPfirst){
    ymat <- t(ymat)
  }




  EigenH5::write_matrix_h5(filename = tempf_X,"dosage",data = xmat,filter="blosc")
  EigenH5::write_matrix_h5(filename = tempf_Y,"trait/ymat",data = ymat,filter="blosc")


  # chunksize <- as.integer(5000)
  # num_chunks <- ceiling(p/chunksize)

  snp_df <- tibble::data_frame(snp_id=1:p) %>%  dplyr::mutate(snp_chunk=ggplot2::cut_number(x=snp_id,n=snp_chunks,labels=F))
  snp_l <-  split(select(snp_df,snp_id,snp_chunk),snp_df$snp_chunk)

  if(SNPfirst){
    snp_lff  <- snp_l %>% purrr::map(~list(subset_rows=.x$snp_id,
                                           filename=tempf_X,
                                           datapath="dosage"))
  }else{
    snp_lff  <- snp_l %>% purrr::map(~list(subset_cols=.x$snp_id,
                                           filename=tempf_X,
                                           datapath="dosage"))
  }

  exp_df <- tibble::data_frame(fgeneid=1:g) %>%  dplyr::mutate(exp_chunk=ggplot2::cut_number(fgeneid,exp_chunks,labels=F))
  exp_l <-  split(select(exp_df,fgeneid,exp_chunk),exp_df$exp_chunk)

  if(EXPfirst){
    exp_lff  <- exp_l %>% purrr::map(~list(subset_rows=.x$fgeneid,
                                           filename=tempf_Y,
                                           datapath="trait/ymat"))
  }else{
    exp_lff  <- exp_l %>% purrr::map(~list(subset_cols=.x$fgeneid,
                                           filename=tempf_Y,
                                           datapath="trait/ymat"))
  }


  snp_exp_lff <-purrr::cross2(snp_lff,exp_lff)

  if(SNPfirst && EXPfirst){
    uh_lff <- map(snp_exp_lff,~list(subset_rows=.x[[1]]$subset_rows,subset_cols=.x[[2]]$subset_rows,filename=uhf,datapath="uh"))
  }
  if((!SNPfirst) && EXPfirst){
    uh_lff <- map(snp_exp_lff,~list(subset_rows=.x[[1]]$subset_cols,subset_cols=.x[[2]]$subset_rows,filename=uhf,datapath="uh"))
  }
  if((SNPfirst) && (!EXPfirst)){
    uh_lff <- map(snp_exp_lff,~list(subset_rows=.x[[1]]$subset_rows,subset_cols=.x[[2]]$subset_cols,filename=uhf,datapath="uh"))
  }
  if((!SNPfirst) && (!EXPfirst)){
    uh_lff <- map(snp_exp_lff,~list(subset_rows=.x[[1]]$subset_cols,subset_cols=.x[[2]]$subset_cols,filename=uhf,datapath="uh"))
  }
  se_lff <-map(uh_lff,~update_list(.,datapath="se"))
  EigenH5::create_matrix_h5(tempf_Z,"uh", numeric(),filter="blosc",dims=as.integer(c(p,g)))
  EigenH5::create_matrix_h5(tempf_Z,"se", numeric(),filter="blosc",dims=as.integer(c(p,g)))
  map_eQTL_chunk_h5(snp_lff,exp_lff,uh_lff,se_lff,cvrt_mat,options=list(SNPfirst=SNPfirst,
                                                                 EXPfirst=EXPfirst))
  cc_uh <- EigenH5::read_matrix_h5(tempf_Z,"uh")
  cc_se <- EigenH5::read_matrix_h5(tempf_Z,"se")
  return(list(uh=cc_uh,se=cc_se))
}



test_that("adjusting for covariates works like matrixeqtl",{

  n <- 10000
  p <- 15
  g <- 6
  ncvrt <- 4
  test_x <- matrix(as.numeric(rbinom(n = n*p,size = 2,prob = 0.1)),p,n)
  test_y <- matrix(as.numeric(rnorm(n = n*g)),g,n)
  bh_se_meqtl <-matrix_eqtl_wrapper(test_x,test_y)
  bh_se_h5 <- eqtl_wrapper(test_x,test_y)
  bh_se_lm <- lm_wrapper(test_x,test_y)
  expect_equal(bh_se_h5,bh_se_meqtl,tolerance=2e-4)
  expect_equal(bh_se_h5,bh_se_lm,tolerance=2e-4)

  cvrt <-matrix(rnorm(n),nrow=n,ncol=1)
  bh_se_meqtl <-matrix_eqtl_wrapper(test_x,test_y,cvrt_mat = t(cvrt))
  bh_se_h5 <- eqtl_wrapper(test_x,test_y,cvrt_mat = cvrt)
  expect_equal(bh_se_h5,bh_se_meqtl,tolerance=2e-4)

  cvrt <-matrix(rnorm(n*ncvrt),nrow=n,ncol=ncvrt)
  bh_se_meqtl <-matrix_eqtl_wrapper(test_x,test_y,cvrt_mat = t(cvrt))
  bh_se_h5 <- eqtl_wrapper(test_x,test_y,cvrt_mat = cvrt)
  expect_equal(bh_se_h5,bh_se_meqtl,tolerance=3e-4)
})


test_that("mapping eQTL works like lm",{
  n <- 10000
  p <- 5
  g <- 3
  test_x <- matrix(as.numeric(rbinom(n = n*p,size = 2,prob = 0.1)),p,n)
  test_y <- matrix(as.numeric(rnorm(n = n*g)),g,n)
  bh_se_h5 <-eqtl_wrapper(test_x,test_y)
  bh_se_lm <- lm_wrapper(test_x,test_y)
  expect_equal(bh_se_h5,bh_se_lm,tolerance=1e-4)
})




test_that("mapping eQTL works like lm with SNPfirst,EXPfirst swapped",{
  n <- 100
  p <- 50
  g <- 30
  test_x <- matrix(as.numeric(rbinom(n = n*p,size = 2,prob = 0.1)),p,n)
  test_y <- matrix(as.numeric(rnorm(n = n*g)),g,n)
  bh_se_lm <- lm_wrapper(test_x,test_y)
  obh_se_h5 <-eqtl_wrapper(test_x,test_y,SNPfirst=TRUE,EXPfirst=TRUE)
  expect_equal(obh_se_h5,bh_se_lm,tolerance=1.5e-2)

  bh_se_h5 <-eqtl_wrapper(test_x,test_y,SNPfirst=TRUE,EXPfirst=FALSE)
  expect_equal(bh_se_h5,obh_se_h5)

  bh_se_h5 <-eqtl_wrapper(test_x,test_y,SNPfirst=FALSE,EXPfirst=TRUE)
  expect_equal(bh_se_h5,obh_se_h5)

  bh_se_h5 <-eqtl_wrapper(test_x,test_y,SNPfirst=FALSE,EXPfirst=FALSE)
  expect_equal(bh_se_h5,obh_se_h5)
})




write_exp_wrapper <-function(xmat,ymat,cvrt_mat=matrix(1,ncol(xmat),1),SNPfirst=TRUE,EXPfirst=TRUE){

  tempf_X <- tempfile()
  tempf_Y <- tempfile()
  tempf_Z <- tempfile()
  p <- nrow(xmat)
  n <- ncol(xmat)
  g <- nrow(ymat)


  stopifnot(ncol(ymat)==n,nrow(cvrt_mat)==n)

  if(!SNPfirst){
    xmat <- t(xmat)
  }
  if(!EXPfirst){
    ymat <- t(ymat)
  }
  EigenH5::write_matrix_h5(filename = tempf_X,"dosage",data = xmat,filter="blosc")
  EigenH5::write_matrix_h5(filename = tempf_Y,"trait/ymat",data = ymat,filter="blosc")


  # chunksize <- as.integer(5000)
  # num_chunks <- ceiling(p/chunksize)

  tibble::data_frame(snp_id=1:p) %>% EigenH5::write_df_h5(tempf_X,"SNPinfo")
  tibble::data_frame(fgeneid=1:g) %>%EigenH5::write_df_h5(tempf_Y,"TraitInfo")
  EigenH5::write_matrix_h5(cvrt_mat,tempf_Z,"covariates")
  return(list(snp_h5=tempf_X,exp_h5=tempf_Y,covar_h5=tempf_Z))

}

test_that("mapping eQTL works like lm with new wrapper",{
  n <- 100
  p <- 50
  g <- 30
  test_x <- matrix(as.numeric(rbinom(n = n*p,size = 2,prob = 0.1)),p,n)
  test_y <- matrix(as.numeric(rnorm(n = n*g)),g,n)
  cvrt_mat=matrix(1,ncol(test_x),1)



  obh_se_h5 <-eqtl_wrapper(test_x,test_y,SNPfirst=TRUE,EXPfirst=TRUE)
  #expect_equal(obh_se_h5,bh_se_lm,tolerance=1.5e-2)

  fl <-write_exp_wrapper(test_x,test_y,SNPfirst=TRUE,EXPfirst=TRUE)
  of <- tempfile()
  map_eqtl_h5(snp_h5 = fl$snp_h5,exp_h5 = fl$exp_h5,covar_h5 = fl$covar_h5,uh_h5 = of)
  res <- list(uh=EigenH5::read_matrix_h5(of,"uh"),se=EigenH5::read_matrix_h5(of,"se"))
  expect_equal(obh_se_h5,res)

  fl <-write_exp_wrapper(test_x,test_y,SNPfirst=TRUE,EXPfirst=FALSE)
  of <- tempfile()
  map_eqtl_h5(snp_h5 = fl$snp_h5,exp_h5 = fl$exp_h5,covar_h5 = fl$covar_h5,uh_h5 = of)
  res <- list(uh=EigenH5::read_matrix_h5(of,"uh"),se=EigenH5::read_matrix_h5(of,"se"))
  expect_equal(obh_se_h5,res)

  fl <-write_exp_wrapper(test_x,test_y,SNPfirst=FALSE,EXPfirst=TRUE)
  of <- tempfile()
  map_eqtl_h5(snp_h5 = fl$snp_h5,exp_h5 = fl$exp_h5,covar_h5 = fl$covar_h5,uh_h5 = of)
  res <- list(uh=EigenH5::read_matrix_h5(of,"uh"),se=EigenH5::read_matrix_h5(of,"se"))
  expect_equal(obh_se_h5,res)

  fl <-write_exp_wrapper(test_x,test_y,SNPfirst=FALSE,EXPfirst=FALSE)
  of <- tempfile()
  map_eqtl_h5(snp_h5 = fl$snp_h5,exp_h5 = fl$exp_h5,covar_h5 = fl$covar_h5,uh_h5 = of)
  res <- list(uh=EigenH5::read_matrix_h5(of,"uh"),se=EigenH5::read_matrix_h5(of,"se"))
  expect_equal(obh_se_h5,res)
})



