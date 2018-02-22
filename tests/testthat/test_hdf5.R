context("hdf5")



test_that("can count up dims in a file",{


  library(purrr)
  n <- 3L
  p <- 4L
  test_mat_1 <-matrix(runif(n*p),n,p)
  test_mat_2 <-matrix(runif(n*p),n,p)

  test_mat_QL <- rerun(11,matrix(runif(p*p),p,p))
  test_mat_uhL <- rerun(11,matrix(runif(n*p),p,n))

out_l <-map2(test_mat_QL,test_mat_uhL,crossprod)

  tfa <-rerun(11,tempfile()) %>% flatten_chr()
  tfb <- rerun(11,tempfile()) %>% flatten_chr()
tfc <- rerun(11,tempfile()) %>% flatten_chr()

iwalk(tfc,~EigenH5::create_matrix_h5(.x,as.character(.y),"result",numeric(),dims=c(p,n)))

pwalk(list(test_mat_QL,tfa,as.character(1:length(test_mat_QL))),~EigenH5::write_matrix_h5(..2,"/",..3,..1))
pwalk(list(test_mat_uhL,tfb,as.character(1:length(test_mat_uhL))),~EigenH5::write_matrix_h5(..2,"/",..3,..1))

  in_df <- tibble::data_frame(filenames=tfa,
                               groupnames="/",
                               datanames=as.character(1:length(test_mat_QL)))
  in2_df <- tibble::data_frame(filenames=tfb,
                              groupnames="/",
                              datanames=as.character(1:length(test_mat_uhL)))
  out_df <-tibble::data_frame(filenames=tfc,
                              groupnames=as.character(1:length(tfc)),
                              datanames="result")
  crossprod_quh_h5(in_df,in2_df,out_df)
  res_mat_l <- imap(tfc,~EigenH5::read_matrix_h5(.x,as.character(.y),"result"))
  expect_equal(res_mat_l,out_l)


  # cum_data_dims(q_dff = sub_df)

})


test_that("crossproduct works",{



    n <- 40L
    p <- 30L
    test_mat_1 <-matrix(runif(n*p),n,p)
    test_mat_2 <-matrix(runif(n*p),n,p)

    tf2 <- tempfile()
    tf3 <- tempfile()
    EigenH5::write_matrix_h5(tf2,"/","test",test_mat_1,doTranspose = F)
    EigenH5::write_matrix_h5(tf3,"/","test",test_mat_2,doTranspose = F)
  tfo <- tempfile()
  Rrres <- crossprod(test_mat_1,test_mat_2)
  crossprod_h5(c(tf2,tf3,tfo),c("/","/","/"),c("test","test","testo"))
  Crres <- EigenH5::read_mat_h5(tfo,"/","testo")
  expect_equal(Rrres,Crres)

  })


})



test_that("matrix (double) roundtrip works with rhdf5",{

  n <- 4L
  p <- 3L
  test_mat <-matrix(1:(n*p)+0.0,n,p)
  tf2 <- tempfile()
  tf3 <- tempfile()
  RcppEigenH5::write_mat_h5(tf2,"/","test",test_mat,doTranspose = T)
  write_mat_h5(tf3,"/","test",test_mat)
  # read_mat <- read_2d_h5(tfile,"/","test",c(0L,0L),c(9L,8L))
  Rread_mat_2 <-RcppEigenH5::read_2d_mat_h5(tf3,"/","test")
  r_read_mat_2 <- read_2d_mat_h5(tf2,"/","test")

  expect_equal(test_mat,Rread_mat_2)
  expect_equal(test_mat,r_read_mat_2)

})


test_that("ld_region splitting works as",{


  num_ld_r <- 20
  region_id <- integer()
  result_matrix <- matrix(0,num_ld_r,3)
  for(i in 1:num_ld_r){
    result_matrix[i,1] <- i
    result_matrix[i,2] <- length(region_id)
    isize <- sample(1:1000,1)
    result_matrix[i,3] <- isize
    region_id <- c(region_id,rep(i,isize))
  }
  N <- 100
  lr <- length(region_id)
  test_mat <- matrix(0.01,lr,N)
  for(i in 1:nrow(result_matrix)){
    row_start <- (result_matrix[i,2]+1)
    row_stop <- row_start+result_matrix[i,3]-1
    row_val <- result_matrix[i,1]
    test_mat[row_start:row_stop,] <- as.numeric(row_val)
  }

  tempf <- tempfile()
  EigenH5::write_vector_h5(filename = tempf,groupname = "SNPinfo",dataname = "region_id",data = region_id)
  write_matrix_h5(filename=tempf,groupname="/",dataname = "dosage",data = test_mat)
  tM <- apply(result_matrix,1,function(x,filename,groupname,dataname,N){
    tX <- read_mat_h5(filename,groupname,dataname,
                      offsets=as.integer(c(x[2],0)),
                      chunksizes = as.integer(c(x[3],N)))
    ttX <- read_ld_chunk_h5(filename,x[1])
    expect_equal(tX,ttX)
  },filename=tempf,groupname="/",dataname="dosage",N=N)


})







test_that("We can read and write dataframes correctly",{


td <- tibble::data_frame(a=1:5,b=letters[1:5],c=runif(5))
tf <- tempfile()
write_df_h5(td,"test",tf)
rtd <- read_df_h5(tf,"test")
expect_equal(td,rtd)
})

test_that("We can read and write vectors correctly",{

  tvec <- runif(30)
  tf <-tempfile()
  write_vec(h5filename = tf,groupname = "test",dataname = "tdata",data = tvec)

  ntvec <- read_vec(tf,"/test/tdata")
  expect_equal(ntvec,tvec)
})



test_that("We can write through groups directly",{

  n <- 4L
  p <- 3L
  test_mat <-matrix(1:(n*p)+0.0,n,p)
  tf3 <- tempfile()
  write_mat_h5(tf3,"testgroup","test",test_mat)
  r_read_mat_2 <- read_2d_mat_h5(tf3,"/testgroup","test")
  expect_equal(test_mat,r_read_mat_2)
  write_mat_h5(tf3,"testgroup/subtest","test",test_mat)
  r_read_mat_2 <- read_2d_mat_h5(tf3,"/testgroup/subtest","test")
  expect_equal(test_mat,r_read_mat_2)
})



test_that("We can read and write integer dataframes",{
  library(dplyr)
  tf <- tempfile()
  tdf <-data_frame(a=1L:10L)
  write_df_h5(tdf,"test_df",tf)
  rdf <- read_df_h5(tf,"test_df")
  expect_equal(tdf,rdf)
  tdf <- mutate(tdf,b=2L:11L)
  tf2 <- tempfile()

  write_df_h5(tdf,"test_df3",tf2)
  rdf <- read_df_h5(tf2,"test_df3")
  expect_equal(tdf,rdf)
})

test_that("We can read and write string dataframes",{
  library(dplyr)
  tf <- tempfile()
  tdf <-data_frame(a=letters)
  write_df_h5(tdf,"test_df",tf)
  rdf <- read_df_h5(tf,"test_df")
  expect_equal(tdf,rdf)
  tdf <- mutate(tdf,b=1:n())
  tf2 <- tempfile()

  write_df_h5(tdf,"test_df3",tf2)
  rdf <- read_df_h5(tf2,"test_df3")
  expect_equal(tdf,rdf)
})

test_that("If I write two strings, one doesn't erase the other",{
  library(dplyr)
  tf <- tempfile()
  tdf <-data_frame(a=letters,b=LETTERS)
  write_df_h5(tdf,"test_df",tf)
  rdf <- read_df_h5(tf,"test_df")
  expect_equal(tdf,rdf)
  tdf <- mutate(tdf,b=1:n())
  tf2 <- tempfile()

  write_df_h5(tdf,"test_df3",tf2)
  rdf <- read_df_h5(tf2,"test_df3")
  expect_equal(tdf,rdf)
})





test_that("Can read specific subcolumns",{
  library(dplyr)
  tf <- tempfile()
  tdf <-data_frame(a=letters,b=LETTERS)
  write_df_h5(tdf,"test_df",tf)
  stdf <- select(tdf,a)
  rdf <- read_df_h5(tf,"test_df",subcols = c("a"))
  expect_equal(stdf,rdf)

  stdf <- select(tdf,b)
  rdf <- read_df_h5(tf,"test_df",subcols = c("b"))
  expect_equal(stdf,rdf)

  stdf <- select(tdf,b,a)
  rdf <- read_df_h5(tf,"test_df",subcols = c("b","a"))
  expect_equal(stdf,rdf)
})




test_that("Can subset dataframes",{
  library(dplyr)
  tf <- tempfile()
  tdf <-data_frame(a=letters,b=LETTERS)
  write_df_h5(tdf,"test_df",tf)
  mfilt <- sample(1:nrow(tdf),4,replace = F)
  stdf <- slice(tdf,mfilt)
  rdf <- read_df_h5(tf,"test_df",filtervec = mfilt)
  expect_equal(stdf,rdf)

})




