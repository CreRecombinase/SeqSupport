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
                            chunksizes=as.integer(c(g,100)))
  EigenH5::write_vector_h5(filename = beta_h5file, groupname = "/", dataname="S",data = runif(p))
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







