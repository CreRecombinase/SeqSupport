context("mapping")

test_that("mapping works like SeqSupport",{

  n <- 100
  p <- 5
  g <- 3
  test_x <- matrix(as.numeric(rbinom(n = n*p,size = 2,prob = 0.1)),p,n)
  test_y <- matrix(as.numeric(rnorm(n = n*g)),g,n)

  tempf_X <- tempfile()
  tempf_Y <- tempfile()
  tempf_Z <- tempfile()
  EigenH5::write_matrix_h5(filename = tempf_X,"/","dosage",data = test_x)
  EigenH5::write_matrix_h5(filename = tempf_Y,"trait","ymat",data = test_y)

  snp_start <- as.integer(seq(0,p-1,by=3))
  snp_stop <- as.integer(pmin(snp_start+3,p))
  snp_size <- as.integer(snp_stop-snp_start)

  snp_dff <- data_frame(filenames=tempf_X,
                        groupnames="/",
                        datanames="dosage",
                       row_offsets=snp_start,
                       row_chunksizes=snp_size,
                       col_offsets=0L,
                       col_chunksizes=as.integer(n))


  gb <- 2
  exp_start <- as.integer(seq(0,g-1,by=gb))
  exp_stop <- as.integer(pmin(exp_start+gb,g))
  exp_size <- as.integer(exp_stop-exp_start)

  exp_dff <- data_frame(filenames=tempf_Y,
                        groupnames="trait",
                        datanames="ymat",
                       row_offsets=exp_start,
                       row_chunksizes=exp_size,
                       col_offsets=0L,
                       col_chunksizes=as.integer(n))

  uh_dff <- cross_df_h5(snp_dff,exp_dff,
                     output_filenames=tempf_Z,
                     output_groupnames="/",
                     output_datanames="uh")
  se_dff <- mutate(uh_dff,datanames="se")

  jj <- 1
  cn <- nrow(snp_dff)*nrow(exp_dff)


  jjm <- matrix(1:cn,nrow(snp_dff),nrow(exp_dff))
  tjj <- 1
  for(i in 1:nrow(exp_dff)){
      for(j in 1:nrow(snp_dff)){
        jj <- jjm[j,i]
        expect_equal(tjj,jj)
        expect_equal(uh_dff$row_offsets[jj],snp_dff$row_offsets[j])
        expect_equal(uh_dff$row_chunksizes[jj],snp_dff$row_chunksizes[j])
        expect_equal(uh_dff$col_offsets[jj],exp_dff$row_offsets[i])
        expect_equal(uh_dff$col_chunksizes[jj],exp_dff$row_chunksizes[i])

        expect_equal(se_dff$row_offsets[jj],snp_dff$row_offsets[j])
        expect_equal(se_dff$row_chunksizes[jj],snp_dff$row_chunksizes[j])
        expect_equal(se_dff$col_offsets[jj],exp_dff$row_offsets[i])
        expect_equal(se_dff$col_chunksizes[jj],exp_dff$row_chunksizes[i])
        tjj <- tjj+1
      }
  }
  EigenH5::create_matrix_h5(tempf_Z,"/","uh", numeric(),dims=as.integer(c(p,g)))
  EigenH5::create_matrix_h5(tempf_Z,"/","se", numeric(),dims=as.integer(c(p,g)))




  map_eQTL_chunk_h5(snp_dff,exp_dff,uh_dff,se_dff)
  cc_uh <- EigenH5::read_matrix_h5(tempf_Z,"/","uh")
  cc_se <- EigenH5::read_matrix_h5(tempf_Z,"/","se")
  cc_bh <- cc_uh*cc_se
  st_x <- apply(t(test_x),2,function(x){x-mean(x)})
  st_y <- apply(t(test_y),2,function(x){x-mean(x)})
  bh <- matrix(0,p,g)
  se <- matrix(0,p,g)
  rbh <- function(x,y){
    sx <-sum(x^2)
  return((t(x)%*%y)/sx)
  }

cpx <- crossprod(st_x,st_y)
ssx <-colSums(st_x^2)
ssy <- colSums(st_y^2)
# rbh <-apply(cpx,2,function(x){x/ssx})
# -2*rbh*cpx+cpx^2/outer(ssx,ssy)
# trbh <-apply(rbh,2,function(x){x^2*(ssx))})
# rse <-sqrt(t(apply(trbh,1,function(x){vy-x})))
for(i in 1:p){
  for(j in 1:g){
    tl <- coef(summary(lm(st_y[,j]~st_x[,i]+0)))
      bh[i,j] <- tl[1]
      se[i,j] <- tl[2]
    }
}
expect_equal(rbh,bh)

  uh <- bh/se
  max(abs(cc_se-se))
  max(abs(cc_uh-uh))
  max(abs(cc_bh-bh))
  expect_equal(cc_se,se)
  expect_equal(cc_uh,uh)
  expect_equal(cc_bh,bh)
})



test_that("crossprod works",{


})

