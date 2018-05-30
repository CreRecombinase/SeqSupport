context("mapping")

test_that("mapping eQTL works like lm",{


  eqtl_wrapper <-function(xmat,ymat){
    tempf_X <- tempfile()
    tempf_Y <- tempfile()
    tempf_Z <- tempfile()
    EigenH5::write_matrix_h5(filename = tempf_X,"/","dosage",data = xmat)
    EigenH5::write_matrix_h5(filename = tempf_Y,"trait","ymat",data = ymat)

    snp_fl <- list(list(filename=tempf_X,datapath="dosage"))
    exp_fl <- list(list(filename=tempf_Y,datapath="trait/ymat"))
    uh_fl <-list(list(filename=tempf_Z,datapath="uh"))
    se_fl <-list(list(filename=tempf_Z,datapath="se"))
    EigenH5::create_matrix_h5(tempf_Z,"/","uh", numeric(),dims=as.integer(c(p,g)))
    EigenH5::create_matrix_h5(tempf_Z,"/","se", numeric(),dims=as.integer(c(p,g)))
    map_eQTL_chunk_h5(snp_fl,exp_fl,uh_fl,se_fl)
    cc_uh <- EigenH5::read_matrix_h5(tempf_Z,"/","uh")
    cc_se <- EigenH5::read_matrix_h5(tempf_Z,"/","se")
    return(list(uh=cc_uh,se=cc_se))
  }
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

  n <- 10000
  p <- 5
  g <- 3
  test_x <- matrix(as.numeric(rbinom(n = n*p,size = 2,prob = 0.1)),p,n)
  test_y <- matrix(as.numeric(rnorm(n = n*g)),g,n)
  bh_se_h5 <-eqtl_wrapper(test_x,test_y)
  bh_se_lm <- lm_wrapper(test_x,test_y)
  expect_equal(bh_se_h5,bh_se_lm,tolerance=1e-4)


})


