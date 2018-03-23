#Code for generating simulations

sim_snp_mat <- function(p=1e4,n=300,af=runif(p)){
  return(matrix(rbinom(n = p*n,size = 2,prob = af),nrow = p,ncol =n))
}


sim_snp_df <- function(p=1e4,break_df=NULL){
data("break_df",package = "LDshrink")
return(tibble::data_frame(chr=sample(1:22,p,replace=T),
                          pos=as.integer(runif(p, 0, .Machine$integer.max)),
                          SNP=paste0("rs",as.integer(runif(p,0,.Machine$integer.max))),
                          allele=replicate(p,paste0(sample(c("A","C","G","T"),2,replace=F),collapse="/")))%>%
  dplyr::distinct() %>% dplyr::arrange(chr,pos)%>% dplyr::mutate(snp_id=1:n()))
}


#' Generate parameters of simulation
#' @param pve true PVE
#' @param bias true value for bias/confounding
#' @param nreps number of repetitions of each combo of parameters
#' @param n sample size
#' @param p number of SNPs
#' @param fgeneid name for each trait (defaults to 1:ntraits)
gen_tparamdf_norm <- function(pve=NULL, bias = 0, nreps, n, p,fgeneid=NULL,sigu=NULL){

  stopifnot(xor(is.null(pve),is.null(sigu)))
  if(!is.null(sigu)){
    pve <-RSSp::calc_pve(sigu = sigu,p_n=p/n)
  }
  stopifnot(all(pve<=1))
  tfgeneid <- fgeneid
  stopifnot(nreps>0)
  if(length(nreps)>1){
    warning("length(nreps)>1, only using first element")
    nreps <- nreps[1]
  }
  tfp <- list(tpve=pve, tbias=bias) %>% purrr::cross_df() %>%
    dplyr::distinct()
  rfp <- dplyr::bind_rows(replicate(nreps, tfp, simplify = F)) %>%
    dplyr::group_by(tpve, tbias) %>%
    dplyr::mutate(replicate=1:n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fgeneid=as.character(1:n()), tsigu=sqrt(n/p*tpve),tbias=tpve*tbias) %>%
    dplyr::select(-replicate)
  if(!is.null(tfgeneid)){
    stopifnot(nrow(rfp)==length(tfgeneid))
    rfp <- dplyr::mutate(rfp,fgeneid=as.character(tfgeneid))
  }
  return(rfp)
}


#' calculate the scale of residuals, given X*Beta and the target PVE
#' @param vy vector with the variance of y (length is equal to the number of traits)
#' @param PVE vector specifying the target PVE (length should be the same as the number of traits)
gen_ti <- function(vy, PVE){
  stopifnot(length(vy)==length(PVE))
  return(vy*(1/PVE-1))
}


#' read log from LD score regression (as implemented in `ldsc`), and parse it so results can be compared to RSSp
#' @param h2lf path (absolute or relative) to `ldsc` (character)
parse_ldsc_h2log <- function(h2lf){
#  library(purrr)
#  library(tidyr)
  h2_dat <- scan(h2lf,what=character(),sep = "\n")
  h2_rowi <- grep("Total Observed scale",h2_dat)
  h2_row <- h2_dat[h2_rowi]
  h2_data <- h2_dat[h2_rowi:length(h2_dat)]
  h2_data <- h2_data[grep("^[^:]+:[^:]+$",h2_data)]
  h2_datd <- purrr::transpose(strsplit(h2_data,split=":"))
  names(h2_datd) <- c("Variable","Value")
  h2_datdf <- tidyr::unnest(as_data_frame(h2_datd)) %>%
    dplyr::mutate(Variable=chartr(" ^","__",Variable),Value=trimws(Value)) %>%
    tidyr::separate(Value,c("Est","SD"),sep = "[ s\\(\\)]+",remove=T,fill="right",convert=T)
  return(h2_datdf)
}


#' Transform simulated standard univariate normal data to multivariate normal, following the RSSp likelihood
#' @param sigu desired sd of the true effect
#' @param bias true value for bias/confounding
#' @param Q matrix of eigenvectors of the LD matrix (must be square, and must have rank equal to nrow(usim))
#' @param D vector of eigenvalues corresponding to Q ( must have length equal to rank of `Q`)
#' @param fgeneid name for the trait (character)
#' @param usim matrix of standard univariate normal data to be transformed, one column for each replicate.
sim_uh_A <- function(sigu,bias,Q,D,fgeneid,usim){
  stopifnot(length(unique(bias))==1,
            length(unique(sigu))==1)
  nd <- sigu^2*D^2+D+bias
  A <- Q %*% (t(Q) * sqrt(pmax(nd, 0)))
  uhmat <- t(usim%*%A)
  colnames(uhmat) <- fgeneid
  return(uhmat)
}

sim_quh <- function(tsigu,tbias,D,fgeneid,p){
  stopifnot(length(tsigu)==length(fgeneid),
            length(tbias)==length(fgeneid))
  return(purrr::map2_dfc(tsigu,tbias,function(sig,a,dvec,p){
    dplyr::data_frame(rnorm(n=p,mean=0,sd=sqrt(sig^2*dvec^2+dvec+a)))
  },dvec=D,p=p) %>% data.matrix() %>% magrittr::set_colnames(as.character(fgeneid)))
}

sim_quh_qu <- function(D,qu,tbias,fgeneid){
  p <- nrow(qu)
  return(purrr::map2_dfc(purrr::array_branch(qu,2),tbias,function(q,a,dvec,p){
    dplyr::data_frame(rnorm(n=p,mean=dvec*q,sd=sqrt(dvec+a)))
  },dvec=D,p=p) %>% data.matrix() %>% magrittr::set_colnames(as.character(fgeneid)))
}



#' "Directly" simulate data from the RSSp likelihood for a single value of `sigu` and `bias`.
#' @param tsigu desired sd of the true effect
#' @param tbias true value for bias/confounding
#' @param Q matrix of eigenvectors of the LD matrix (must be square, and must have rank equal to nrow(usim))
#' @param D vector of eigenvalues corresponding to Q ( must have length equal to rank of `Q`)
#' @param fgeneid name for the trait (character)
#' @param snp_id vector of IDs for the SNPs
sim_uh_quh_dir_u <- function(tsigu,tbias,fgeneid=NULL,Q,D,seed=NULL,snp_id,umat){
  stopifnot(length(tsigu)==length(tbias),
            length(tbias)==length(fgeneid))
  p <- length(snp_id)
  fgeneid <- as.character(fgeneid)
  qu <-crossprod(Q,umat)
  quh <-sim_quh_qu(D,qu,tbias,fgeneid)
  uh <-Q%*%quh
  colnames(uh) <- fgeneid
  return(list(quh=quh,snp_id=snp_id,uh=uh))
}





#' "Directly" simulate data from the RSSp likelihood for a single value of `sigu` and `bias`.
#' @param tsigu desired sd of the true effect
#' @param tbias true value for bias/confounding
#' @param Q matrix of eigenvectors of the LD matrix (must be square, and must have rank equal to nrow(usim))
#' @param D vector of eigenvalues corresponding to Q ( must have length equal to rank of `Q`)
#' @param fgeneid name for the trait (character)
#' @param snp_id vector of IDs for the SNPs
sim_uh_quh_dir <- function(tsigu,tbias,fgeneid=NULL,Q,D,seed=NULL,snp_id){
stopifnot(length(tsigu)==length(tbias),
          length(tbias)==length(fgeneid))

  p <- length(snp_id)
  fgeneid <- as.character(fgeneid)
  quh <-sim_quh(tsigu,tbias,D,fgeneid,p)
  uh <-Q%*%quh
  colnames(uh) <- fgeneid
  return(list(quh=quh,snp_id=snp_id,uh=uh))
}



#' "Directly" simulate data from the RSSp likelihood for a range of values of sigu and bias
#' @param tsigu desired sd of the true effect (vector of length `g``, where `g` is the number of traits to be simulated)
#' @param tbias true value for bias/confounding
#' @param Q matrix of eigenvectors of the LD matrix (must be square, and must have rank equal to nrow(usim))
#' @param D vector of eigenvalues corresponding to Q ( must have length equal to rank of `Q`)
#' @param fgeneid name for the trait (character)
#' @param snp_id vector of IDs for the SNPs
sim_quh_dir_df <- function(tparam_df,D,seed=NULL,snp_id){
  tparam_dfl <- split(tparam_df, tparam_df$fgeneid)
  quh <- do.call(cbind,lapply(tparam_dfl,function(l,D,seed,snp_id){
    sim_quh_dir(
      tsigu=l$tsigu,
      tbias=l$tbias,
      fgeneid=l$fgeneid,
      D=D,
      seed=seed,
      snp_id=snp_id)[["quh"]]
  },D=D,seed=seed,snp_id=snp_id))
  return(list(quh=quh,D=D,snp_id=snp_id))
}





#' Estimate parameters from simulated data, saved in the specified format
#' @param resl list with the following elements
#' `quh_mat`, a `p`x`g` matrix (where `p` is number of SNPs and `g` is number of traits) obtained by an
#' operation equivalent to to `crossprod(Q,uh)`
#' `tparam_df` a dataframe (of the type generated by the function `gen_tparamdf_norm` with true parameter values (`NA`s can also be used in the case of real data)
#' @param Ql A list of matrices, each element representing the eigenvectors of a square block in the block diagonal approximation of the LD matrix (blocks need not be of equal size)
#' @param D A vector of the eigenvectors of the LD matrix (must be length `p`)
#' @param doConfound a logical indicating whether to use the two parameter model (`doConfound=T`) or the one parameter model without confounding
#' @param log_params a logical indicating whether to optimize parameters in log space or not (not recommended)
#' @param a_bounds a length two vector indicating the bounds of the parameter `a` to use when optimizing, `a` is a measure of confounding
#' @param sigu_bounds a length two vector indicating bounds in the parameter `sigu` to use when optimizing `sigu` is directly related to PVE
est_sim <- function(resl,Ql=NULL,D=NULL,doConfound=T,log_params=F,useGradient=T,bias_bounds=c(0,.3),pve_bounds=c(1e-4,1)){

  stopifnot(!is.null(D))
  if(is.null(Ql)){
    stopifnot(!is.null(resl$quh_mat))
  }
  if(is.null(resl$quh_mat)){
    resl$quh_mat <- RSSp::quh_mat(Ql,resl$bias_uh_mat)
  }


  rss_res <- purrr::cross(list(doConfound=doConfound,log_params=log_params,useGradient=useGradient)) %>%
    purrr::invoke_map_dfr(
      RSSp::RSSp_run_mat_quh,
      .,
      quh_mat=resl$quh_mat,
      D=D,
      n=resl$n,
      a_bounds=bias_bounds,
      pve_bounds=pve_bounds) %>% dplyr::inner_join(resl$tparam_df,by="fgeneid")
  return(rss_res)
}

RSSp_oracle <-function(quh,dvec,tsigu,tbias,fgeneid=NULL){

  return(tibble::data_frame(fgeneid=fgeneid,
                     oracle_lnZ=-RSSp::evd_dnorm(par=c(tsigu^2,tbias),
                                           dvec=dvec,
                                           quh=quh),
                     tsigu_cov=-1/RSSp::RSSp_hess(c(tsigu,tbias),dvec=dvec,quh=quh)[1,1]))
}




sim_S <- function(index,x,sigu,outfile_h5=NULL){
  #I'm adding a sneaky bit in here to deal with the (extremely) rare
  #event where all individuals are hets at a particular locus
  #what I'm doing is simply assigning the variant a beta of 0
  n <- nrow(x)
  p <- ncol(x)
  sx <- scale(x,center=T,scale=F)
  u_mat <- sapply(sigu,function(tsigu,p){rnorm(n=p,mean=0,sd=tsigu)},p=p)
  S <- 1/sqrt(n)*1/apply(sx,2,sd)
  beta <- u_mat*S
  beta[!is.finite(beta)]<-0
  ty <- sx%*%beta
  if(!is.null(outfile_h5))
  {
    tbeta <- t(beta)
    S <- c(S)
    h5write(obj=tbeta,
            file=outfile_h5,
            name="beta",
            start=c(1,index))
    h5write(obj=S,
            file=outfile_h5,
            name="S",
            start=index)
  }
  return(ty)
}


#' Generate simulated values of beta, given PVE etc.
#' @param gds an open SeqArray gds file handle
#' @param tparam_df A dataframe with the simulation parameter values (see `gen_tparamdf_norm`)
#' @param seed a random seed
#' @param chunksize number of values to create at a time
sim_beta_gds <- function(gds,tparam_df,seed=NULL,chunksize=10000){

  sim_beta <- function(x,sigu,fgeneid){
    rsid=x$rsid
    allele=x$allele
    chrom=x$chrom
    pos=x$pos
    x_g <- x$x

    sx <- scale(x_g,center=T,scale=F)
    #I'm adding a sneaky bit in here to deal with the (extremely) rare
    #event where all individuals are hets at a particular locus
    #what I'm doing is simply assigning the variant a beta of 0
    n <- nrow(sx)
    p <- ncol(sx)
    u_mat <- sapply(sigu,function(tsigu,p){rnorm(n=p,mean=0,sd=tsigu)},p=p)
    S <- 1/sqrt(n)*1/apply(sx,2,sd)
    beta <- u_mat*S
    beta[!is.finite(beta)]<-0
    colnames(beta) <- fgeneid
    tdf <- tibble::data_frame(SNP=rsid,
                              allele=allele) %>%
      separate(allele,into = c("A1","A2"),sep = ",") %>%
      mutate(snp_id=1:n())

    tdf<- as_data_frame(beta) %>% mutate(snp_id=1:n()) %>%
      gather(fgeneid,beta,-snp_id) %>% inner_join(tdf) %>%
      select(SNP,A1,A2,beta,fgeneid)
    return(tdf)
  }


  result <- seqBlockApply(gds,c(x="$dosage",rsid="annotation/id",
                                allele="allele",
                                chrom="chromosome",
                                pos="position"),
                          FUN=sim_beta,sigu=tparam_df$tsigu,
                          fgeneid=tparam_df$fgeneid,
                          margin="by.variant",
                          as.is = "list",
                          .progress = T,bsize = chunksize,parallel=F)
  res <- bind_rows(result)
  return(res)

}

#' Find number of SNPs in the dataset
#' @param gds A SeqArray gds object
calc_p <- function(gds){
  length(SeqArray::seqGetData(gds,"variant.id"))
}

calc_p_h5 <- function(h5){
EigenH5::get_dims_h5(h5[1],h5[2],h5[3])[1]
}
calc_N_h5 <- function(h5){
  EigenH5::get_dims_h5(h5[1],h5[2],h5[3])[2]
}

#' Find number of individuals in the dataset
#' @param gds A SeqArray gds object
calc_N <- function(gds){
  length(SeqArray::seqGetData(gds,"sample.id"))
}






gen_ty_h5 <- function(snp_df,snp_h5file,beta_h5file,tparam_df,chunksize=1000){

  p <- as.integer(nrow(snp_df))
  N <-as.integer(unique(tparam_df$n))
  stopifnot(all(tparam_df$p==p))
  g <- as.integer(nrow(tparam_df))
  stopifnot(length(N)==1)
  dims_beta <-as.integer(c(g,p))
  num_chunks <- ceiling(p/chunksize)
  input_dff <- dplyr::mutate(snp_df,nchunk_id=sort(as.integer(gl(n = num_chunks,k=chunksize,length = p)))) %>%
    EigenH5::split_chunk_df(pos_id=snp_id,group_id=nchunk_id) %>%
    dplyr::mutate(chunk_group=nchunk_id) %>%
    dplyr::filter(!is.na(chunk_group))
  stopifnot(sum(input_dff$row_chunksizes)==p)
  d_dims <-EigenH5::get_dims_h5(snp_h5file,"/","dosage")
  SNPfirst <- d_dims[2]==N
  if(SNPfirst){
    mr <- filter(input_dff,row_offsets+row_chunksizes>d_dims[1])
  }else{
    mr <- filter(input_dff,row_offsets+row_chunksizes>d_dims[2])
  }
  stopifnot(nrow(mr)==0)

  stopifnot(nrow(input_dff)==num_chunks)
  beta_dff <- data_frame(id=1:p,nchunk_id=sort(as.integer(gl(n = num_chunks,k=chunksize,length = p)))) %>%
    EigenH5::split_chunk_df(pos_id=id,group_id=nchunk_id) %>%
    dplyr::mutate(chunk_group=nchunk_id,filenames=beta_h5file,groupnames="/",datanames="Beta",row_offsets=0L,row_chunksizes=g)
  stopifnot(nrow(beta_dff)==nrow(input_dff))
  stopifnot(sum(beta_dff$col_chunksizes)==p)

  data_dff <- dplyr::mutate(input_dff,filenames=snp_h5file,groupnames="/",datanames="dosage")
  if(SNPfirst){
    data_dff <- dplyr::mutate(data_dff,col_offsets=0L,col_chunksizes=N)
  }else{
    data_dff <- dplyr::mutate(data_dff,row_offsets=0L,row_chunksizes=N)
  }
  S_dff <- dplyr::mutate(beta_dff,filenames=beta_h5file,groupnames="/",datanames="S",
                         row_offsets=col_offsets,
                         row_chunksizes=col_chunksizes,
                         col_offsets=0L,col_chunksizes=1L)
  out_dff <- dplyr::bind_rows(beta_dff,S_dff)
  EigenH5::create_matrix_h5(filename = beta_h5file,
                            groupname = "/",
                            dataname = "Beta",
                            data = numeric(),
                            dims = dims_beta,
                            doTranspose = F,
                            chunksizes=as.integer(c(g,100)))
  EigenH5::create_vector_h5(filename = beta_h5file, groupname = "/", dataname="S",dimension=as.integer(p),chunksize=as.integer(100))
  ymat <- simulate_y_h5(data_dff,out_dff,p,N,g,tparam_df$tsigu)
  return(ymat)
}


gen_sim_resid <- function(ty,tparam_df){
  vy <- apply(ty, 2, var)
  n <- nrow(ty)
  stopifnot(ncol(ty)==nrow(tparam_df))
  residvec <- gen_ti(vy, tparam_df$tpve)
  residmat <- sapply(residvec, function(ti, n){rnorm(n=n, mean=0, sd=sqrt(ti))}, n=n)
  ymat <- scale(ty+residmat, center=T, scale=F)
  return(ymat)
}


gen_sim_phenotype_h5 <- function(snp_df,snp_h5file,beta_h5file,tparam_df,chunksize=1000){
  # if(!is.null(seed)){
  #   set.seed(seed)
  # }
  # 1000
  ty <- gen_ty_h5(snp_df,snp_h5file,beta_h5file,tparam_df,chunksize=chunksize)
  # h5file=h5file,tparam_df=tparam_df,seed=seed,chunksize=chunksize,betamat=betamat,beta_h5file=beta_h5file)

  return(gen_sim_resid(ty,tparam_df))
}




read_SNPinfo_ldsc_gwas <- function(gds,zmat,N=NULL){
  if(is.null(N)){
    N <- calc_N(gds)

    if(LDshrink::is_haplo(gds)){
      N <- N/2
    }
  }
  tdf <- tibble::data_frame(SNP=seqGetData(gds,var.name="annotation/id"),
                            allele=seqGetData(gds,var.name="allele")) %>%
    tidyr::separate(allele,into = c("A1","A2"),sep = ",") %>%
    dplyr::mutate(N=N,snp_id=1:n())
  tdf<- dplyr::as_data_frame(zmat) %>% dplyr::mutate(snp_id=1:n()) %>%
    tidyr::gather(fgeneid,Z,-snp_id) %>%
    dplyr::inner_join(tdf) %>%
    dplyr::select(SNP,A1,A2,Z,N,fgeneid)
  return(split(select(tdf,-fgeneid),tdf$fgeneid))
}







