context("simulations")

test_that("pve and sigu are calculated consistently",{
  p <- 100
  n <- 10
  tpve <- 0.4
  sigu <- RSSp::calc_sigu(tpve,p/n)
  pve <- RSSp::calc_pve(sigu,p/n)
  expect_equal(pve,tpve)

  tsigu <- .3
  pve <- RSSp::calc_pve(tsigu,p/n)
  sigu <-RSSp::calc_sigu(pve,p/n)
  expect_equal(tsigu,sigu)
})

test_that("normal simulations are kind of normal",{
p <- 100000
tsigu <-c(1,1.2,1.5,2,3)
U <- SeqSupport::sim_U(n = p,tsigu)$U
expect_equal(dim(U),c(p,length(tsigu)))
expect_true((!is.na(mean(U))))
expect_gt(var(U[,2]),var(U[,1]))
expect_equal(apply(U,2,sd),tsigu,tolerance = 0.1)
expect_equal(colMeans(U,2),rep(0,length(tsigu)),tolerance=0.1)
})

test_that("normal simulations (using tbb) are kind of normal",{
  p <- 10000
  tsigu <-c(1,1.2,1.5,2,3)


  U <-  SeqSupport::sim_U2(n = p,tsigu)

  #                U= SeqSupport::sim_U(n = p,tsigu))


  expect_equal(dim(U),c(p,length(tsigu)))
  expect_true((!is.na(mean(U))))
  expect_gt(var(U[,2]),var(U[,1]))
  expect_equal(apply(U,2,sd),tsigu,tolerance = 0.1)
  expect_equal(colMeans(U,2),rep(0,length(tsigu)),tolerance=0.1)
})

test_that("pve and sigu are calculated consistently in gen_tparam_df",{


  p <- 100
  n <- 10
  tpve <- 0.4
  tp_df <- gen_tparamdf_norm(pve=tpve,bias=0,nreps=1,n=n,p=p)
  expect_equal(tp_df$tpve,tpve)
  expect_equal(tp_df$tsigu,RSSp::calc_sigu(tpve,p/n))

})


test_that("AF is calculated correctly",{
  n <- 10000
  p <- 5
  g <- 3
  tp <- runif(p)
  tX <- sapply(tp,rbinom,n=n,size=2)
  tempf <- tempfile()
  EigenH5::write_matrix_h5(filename = tempf,groupname = "/",dataname = "dosage",data=t(tX))
  fl <-list(list(filename=tempf,datapath="dosage",subset_rows=1:3),
            list(filename=tempf,datapath="dosage",subset_rows=4:5))
  res_af <- calc_af_h5(fl)
  expect_equal(res_af,alt_AF(tX))
  expect_equal(res_af,tp,tolerance=0.01)
  # test_x <- matrix(as.numeric(rbinom(n = n*p,size = 2,prob = 0.1)),p,n)
})

# test_that("Can simulate samples from EVD",{
#
#   n <- 10000
#   p <- 3
#
#   tX <- scale(matrix(rnorm(n*p), n, p), center = T, scale = F)
#   tS <- cov(tX)
#   tC <- eigen(tS)
#   nv <- matrix(rnorm(n*p),p,n)
#   # r_sim <- rnorm_evd_int(tC$vectors, sqrt(tC$values), nv)
#   # c_sim <- evd_rnorm_i(tC$vectors, sqrt(tC$values), nv)
#   cS <- cov(t(c_sim))
#   rS <- cov(t(r_sim))
#   expect_equal(r_sim,c_sim)
#   expect_lt(max(abs(cS - tS))/max(abs(tS)), 0.1)
# })


