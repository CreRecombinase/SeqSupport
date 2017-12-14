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


test_that("svd",{
  xf <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/h5/sub_seq_train_hapmap_geno.h5"
      yf <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/RSSp_sim_gwas_pheno/sub_NoConfoundSmall_trait.h5"

    xf <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/h5/scz2_seq_train_hapmap_geno.h5"
    yf <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/RSSp_sim_gwas_pheno/scz2_NoConfound_trait.h5"

  X <- RcppEigenH5::read_2d_mat_h5(xf,"/","dosage")
  y <- RcppEigenH5::read_2d_mat_h5(yf,"trait","ymat")
  si <- read_df_h5(yf,"SimulationInfo")
stopifnot(nrow(y)==ncol(X))
  # We want to scale the ROWS not the columns (because it's a kinship matrix NOT an LD matrix)
  p <- nrow(X)
  N <- ncol(X)
  X <- t(X)
  X <- scale(X,center=T,scale=F)
  K <- tcrossprod(X)/p
  svdX <- svd(X/(sqrt(p)))
  D <- svdX$d
  U <- svdX$u
  V <- svdX$v
  rm(svdX)
  gc()
  sa <- runif(1)
  H <- diag(N) + sa*K
  sH <- solve(H,y[,1])
  tsH <- solve(H)%*%y[,1]
  expect_equal(c(tsH),c(sH))
  expect_equal(c(H%*%sH),y[,1])




  svd_sH <- U%*%diag(1/(sa*D^2+1))%*%(t(U)%*%y[,1])
  expect_equal(t(U)%*%diag(1/(D^2+1))%*%U,solve(K+diag(N)))
  expect_equal(t(U)%*%diag(1/(sa*D^2+1))%*%U,solve(H))
  expect_equal(sum(svd_sH),sum(tsH))
  expect_equal(U%*%diag(1/(D^2))%*%t(U),s_K)
  expect_equal(U%*%diag(D^2)%*%t(U),K)
  expect_equal(Q%*%diag(eD)%*%t(Q),K)
  expect_equal(sH,c(svd_sH))
  expect_equal(sH,c(evd_sH))
  svd_b <- sum(y[,1]*U%*%diag(1/(sa*D^2+1))%*%(t(U)%*%y[,1]))
  b        <- sum(y[,1]*solve(H,y[,1]))
  expect_equal(svd_b,b)
  td <- determinant(b*H,logarithm = T)$modulus
  otd <- determinant(H,logarithm = T)
  expect_equivalent(,td)
  expect_equal(log(b)*N+otd$modulus,td)
  p_det <- as.numeric(-determinant(b*H,logarithm = TRUE)$modulus/2)
  np_det <- -(N*log(b)+determinant(H,logarithm = T)$modulus)/2
  expect_equal(np_det,p_det)
  svd_det <- -(N*log(b)+sum(log(sa*D^2+1)))/2
  expect_equal(svd_det,p_det)
  expect_equal(svd_b,b)
  expect_equivalent(U%*%diag(D)%*%t(V),as.matrix(X/sqrt(p)))
  expect_equivalent((U)%*%diag(D^2)%*%t(U),K)
  sx <- sum(apply(X,2,sd)^2)
  reg_mod <- polygenic.model(X,y[,1],0.4)
  p <- nrow(V)

  svd_mod <- polygenic.model_svd(U = U,D = D,y = y[,1],h = 0.4,sx = sum(apply(X,2,sd)^2),p = p)
  expect_equal(reg_mod,svd_mod)

all_res <- apply(y,2,function(yv,x){
    optimise(optim_pg,interval = c(1e-4,1),X=X,y=yv)

  },x=X)
all_r <- map_df(all_res,as_data_frame) %>% mutate(fgeneid=as.character(1:n()))
inner_join(all_r,si) %>% ggplot(aes(x=tpve,y=minimum))+geom_point()+geom_smooth(method="lm")+geom_abline(slope = 1,intercept = 0)
})


