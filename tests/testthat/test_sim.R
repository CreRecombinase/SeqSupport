context("simulations")

test_that("pve and sigu are calculated consistently",{


  p <- 100
  n <- 10
  tpve <- 0.4
  sigu <- calc_sigu(tpve,p/n)
  pve <- calc_pve(sigu,p/n)
  expect_equal(pve,tpve)

  tsigu <- .3
  pve <- calc_pve(tsigu,p/n)
  sigu <-calc_sigu(pve,p/n)
  expect_equal(tsigu,sigu)
})



test_that("pve and sigu are calculated consistently in gen_tparam_df",{


  p <- 100
  n <- 10
  tpve <- 0.4
  tp_df <- gen_tparamdf_norm(pve=tpve,bias=0,nreps=1,n=n,p=p)
  expect_equal(tp_df$tpve,tpve)
  expect_equal(tp_df$tsigu,calc_sigu(tpve,p/n))

})
