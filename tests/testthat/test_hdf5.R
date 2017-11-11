context("hdf5")





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




