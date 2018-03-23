context("Matching functions")

test_that("matching SNPs to regions works with one chromosome",{
  library(tidyverse)
  nregs <- 10
  reg_test <- data_frame(chr = 3L,start = as.integer(seq(1,1e6,length.out = 10))) %>%
    mutate(stop = start + sample(10:1000,size = n(),replace = T),
           range_id = paste0("range: ",1:n()),
           num_snps = sample(0:5,n(),replace = T))
  good_snps <- split(reg_test,reg_test$range_id) %>% map(function(x){
    tret <- bind_rows(replicate(x$num_snps,data_frame(chr=x$chr,pos=as.integer(runif(n=1,x$start,x$stop))),simplify=F))
    if(nrow(tret)>0){
      return(tret)
    }else{
      return(NULL)
    }
  }) %>% bind_rows() %>% mutate(snp_id = paste0("SNP: ",1:n()))
  test_map <- match_SNP(reg_test,good_snps)
  snp_ct <- group_by(test_map,range_id) %>% summarise(n_snps=n()) %>% inner_join(distinct(test_map,range_id,num_snps))
  expect_equal(snp_ct$num_snps,snp_ct$n_snps)
})

test_that("boost_transpose",{
am <- matrix(as.numeric(1:18),9,2)
oam <-matrix(as.numeric(1:18),9,2)
tam <- matrix(as.numeric(1:18),9,2,byrow = T)
test_transpose_matrix(am)

})


test_that("convert to and from GDS",{

  snp_df <- sim_snp_df()
  X <- sim_snp_mat()
  tf <- snpgdsR2SNP(X,snp_df)
  ttf <- tempfile()
  tttf <- tempfile()
  SeqArray::seqSNP2GDS(tf,ttf)
  tgds <- SeqArray::seqOpen(ttf)
  gds2hdf5(gdsfile = tgds,hdf5file = tttf)
  n_snp_df <- EigenH5::read_df_h5(tttf,"SNPinfo")
})
