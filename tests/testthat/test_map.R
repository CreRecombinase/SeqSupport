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

eqtl_wrapper <-function(xmat,ymat,cvrt_mat=matrix(1,ncol(xmat),1)){
  tempf_X <- tempfile()
  tempf_Y <- tempfile()
  tempf_Z <- tempfile()
  EigenH5::write_matrix_h5(filename = tempf_X,"/","dosage",data = xmat,filter="blosc")
  EigenH5::write_matrix_h5(filename = tempf_Y,"trait","ymat",data = ymat,filter="blosc")

  snp_fl <- list(list(filename=tempf_X,datapath="dosage"))
  exp_fl <- list(list(filename=tempf_Y,datapath="trait/ymat"))
  uh_fl <-list(list(filename=tempf_Z,datapath="uh"))
  se_fl <-list(list(filename=tempf_Z,datapath="se"))
  # EigenH5::start_blosc()
  EigenH5::create_matrix_h5(tempf_Z,"/","uh", numeric(),filter="blosc",dims=as.integer(c(p,g)))
  EigenH5::create_matrix_h5(tempf_Z,"/","se", numeric(),filter="blosc",dims=as.integer(c(p,g)))
  map_eQTL_chunk_h5(snp_fl,exp_fl,uh_fl,se_fl,cvrt_mat)
  cc_uh <- EigenH5::read_matrix_h5(tempf_Z,"/","uh")
  cc_se <- EigenH5::read_matrix_h5(tempf_Z,"/","se")
  return(list(uh=cc_uh,se=cc_se))
}



test_that("adjusting for covariates works like matrixeqtl",{

  n <- 10000
  p <- 5
  g <- 3
  ncvrt <- 4
  test_x <- matrix(as.numeric(rbinom(n = n*p,size = 2,prob = 0.1)),p,n)
  test_y <- matrix(as.numeric(rnorm(n = n*g)),g,n)
  bh_se_meqtl <-matrix_eqtl_wrapper(test_x,test_y)
  bh_se_h5 <- eqtl_wrapper(test_x,test_y)
  bh_se_lm <- lm_wrapper(test_x,test_y)
  expect_equal(bh_se_h5,bh_se_meqtl,tolerance=2e-4)

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


