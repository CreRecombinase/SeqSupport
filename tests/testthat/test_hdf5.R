context("hdf5")


#
# test_that("can count up dims in a file",{
#
#
#   library(purrr)
#   n <- 3L
#   p <- 4L
#   test_mat_1 <-matrix(runif(n*p),n,p)
#   test_mat_2 <-matrix(runif(n*p),n,p)
#
#   test_mat_QL <- rerun(11,matrix(runif(p*p),p,p))
#   test_mat_uhL <- rerun(11,matrix(runif(n*p),p,n))
#
# out_l <-map2(test_mat_QL,test_mat_uhL,crossprod)
#
#   tfa <-rerun(11,tempfile()) %>% flatten_chr()
#   tfb <- rerun(11,tempfile()) %>% flatten_chr()
# tfc <- rerun(11,tempfile()) %>% flatten_chr()
#
# iwalk(tfc,~EigenH5::create_matrix_h5(.x,as.character(.y),"result",numeric(),dims=c(p,n)))
#
# pwalk(list(test_mat_QL,tfa,as.character(1:length(test_mat_QL))),~EigenH5::write_matrix_h5(..2,"/",..3,..1))
# pwalk(list(test_mat_uhL,tfb,as.character(1:length(test_mat_uhL))),~EigenH5::write_matrix_h5(..2,"/",..3,..1))
#
#   in_df <- tibble::data_frame(filenames=tfa,
#                                groupnames="/",
#                                datanames=as.character(1:length(test_mat_QL)))
#   in2_df <- tibble::data_frame(filenames=tfb,
#                               groupnames="/",
#                               datanames=as.character(1:length(test_mat_uhL)))
#   out_df <-tibble::data_frame(filenames=tfc,
#                               groupnames=as.character(1:length(tfc)),
#                               datanames="result")
#   crossprod_quh_h5(in_df,in2_df,out_df)
#   res_mat_l <- imap(tfc,~EigenH5::read_matrix_h5(.x,as.character(.y),"result"))
#   expect_equal(res_mat_l,out_l)
#
#
#   # cum_data_dims(q_dff = sub_df)
#
# })
#
#
#
#
#
#



