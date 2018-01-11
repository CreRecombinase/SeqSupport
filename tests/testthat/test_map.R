context("mapping")

test_that("mapping works like SeqSupport",{



bdf <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/h5/bd_seq_hapmap_geno.h5"
tX <- EigenH5::read_mat_h5(bdf,"/","dosage")
  n <- 200
  p <- 5
  g <- 3

  test_X <- matrix(as.numeric(rbinom(n = n*p,size = 2,prob = 0.1)),n,p)
  test_y <- matrix(as.numeric(rnorm(n = n*g)),n,g)
  tempf_X <- tempfile()
  tempf_Y <- tempfile()
  tempf_Z <- tempfile()
  EigenH5::write_mat_h5(filename = tempf_X,"/","dosage",data = t(test_X))
  EigenH5::write_mat_h5(filename = tempf_Y,"/","ymat",data = t(test_y))

  map_eQTL2_h5(c(tempf_X,"/","dosage"),
              c(tempf_Y,"/","ymat"),
              out_path = c(tempf_Z,"eQTL"),
              SNP_chunksize = -1,EXP_chunksize = -1)
  EigenH5::map_eQTL_h5(c(tempf_X,"/","dosage"),
                       c(tempf_Y,"/","ymat"),
                       out_path = c(tempf_Z,"eQTL"),
                       SNP_chunksize = -1,EXP_chunksize = -1)
  c_uh <- read_mat_h5(tempf_Z,"eQTL","uh")
  c_se <- read_mat_h5(tempf_Z,"eQTL","se")
  c_bh <- c_uh*c_se
  st_X <- apply(test_X,2,function(x){x-mean(x)})
  st_y <- apply(test_y,2,function(x){x-mean(x)})
  bh <- matrix(0,p,g)
  se <- matrix(0,p,g)
  for(i in 1:p){
    for(j in 1:g){
      tl <- coef(summary(lm(st_y[,j]~st_X[,i]+0)))
      bh[i,j] <- tl[1]
      se[i,j] <- tl[2]
    }
  }
  uh <- bh/se
  max(abs(c_se-se))
  max(abs(c_uh-uh))
  max(abs(c_bh-bh))
  expect_equal(c_se,se,tolerance=2.9e-3)
  expect_equal(c_uh,uh,tolerance=0.01)
  expect_equal(c_bh,bh)
})
