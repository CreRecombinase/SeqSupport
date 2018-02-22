context("mapping")

test_that("mapping works like SeqSupport",{

  n <- 100
  p <- 5
  g <- 3

  test_x <- matrix(as.numeric(rbinom(n = n*p,size = 2,prob = 0.1)),n,p)
  test_y <- matrix(as.numeric(rnorm(n = n*g)),n,g)

  tempf_X <- tempfile()
  tempf_Y <- tempfile()
  tempf_Z <- tempfile()
  EigenH5::write_matrix_h5(filename = tempf_X,"/","dosage",data = t(test_x))
  # reg_v <- as.integer(rep(1,p))
  # EigenH5::write_vector_h5(tempf_X,"SNPinfo","region_id",reg_v)
  EigenH5::write_matrix_h5(filename = tempf_Y,"trait","ymat",data = test_y)
  y_df <- chunk_df_h5(filename=tempf_Y,groupname="trait",dataname="ymat")

  X_df <- chunk_df_h5(filename=tempf_X,groupname="/",dataname="dosage",chunksize_row = 3)
  out_dff <-gen_outer_df(y_df,X_df)


  map_eQTL_chunk_h5(tempf_X,tempf_Y,tempf_Z)
  cc_uh <- EigenH5::read_matrix_h5(tempf_Z,"1","uh")
  cc_se <- EigenH5::read_matrix_h5(tempf_Z,"1","se")
  cc_bh <- cc_uh*cc_se
  st_x <- apply(test_x,2,function(x){x-mean(x)})
  st_y <- apply(test_y,2,function(x){x-mean(x)})
  bh <- matrix(0,p,g)
  se <- matrix(0,p,g)

  for(i in 1:p){
    for(j in 1:g){
      tl <- coef(summary(lm(st_y[,j]~st_x[,i]+0)))
      bh[i,j] <- tl[1]
      se[i,j] <- tl[2]
    }
  }
  uh <- bh/se
  max(abs(cc_se-se))
  max(abs(cc_uh-uh))
  max(abs(cc_bh-bh))
  expect_equal(cc_se,se)
  expect_equal(cc_uh,uh)
  expect_equal(cc_bh,bh)
})


test_that("chunk crossprod works",{


  p <- 50
  g <- 3
  chunks <- 50
  bfmat <- matrix(0,chunks*p,g)
  il <- split(1:(chunks*p),cut(1:chunks*p,breaks = chunks,labels = F))
  tfa <- tempfile()
  tfb <- tempfile()
Ql <- list()
uhl <- list()
quhl <- list()
  for(i in 1:chunks){
    Ql[[i]] <-matrix(runif(p*p),p,p)
    uhl[[i]] <- matrix(runif(p*g),p,g)
    quhl[[i]] <- crossprod(Ql[[i]],uhl[[i]])
    EigenH5::write_matrix_h5(filename = tfa,groupname = paste0("EVD/",as.character(i)),dataname = "Q",data = Ql[[i]])
    EigenH5::write_matrix_h5(filename = tfb,groupname = as.character(i),dataname = "uh",data = uhl[[i]])
  }
bfmat <- data.matrix(map_df(quhl,as_data_frame))
attr(bfmat,"dimnames") <- NULL
  tfo <- tempfile()
  crossprod_chunk_h5(infile_1 = tfa,
                     infile_2 = tfb,in_groups = as.character(1:50),
                     outfile = tfo,indata_1 = "Q",indata_2 = "uh",outgroup = "/",outdata = "quh")
  rquh <- EigenH5::read_matrix_h5(tfo,"/","quh")
  expect_equal(bfmat,rquh)


})
