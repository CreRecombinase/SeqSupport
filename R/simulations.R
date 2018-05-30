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
    retvec <- vy*(1/PVE-1)
    retvec[PVE==0] <- 1
    return(retvec)
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


calc_S <- function(centered_X,margin=2){
  if(margin==2){
    N <- nrow(centered_X)
  }else{
    N <- ncol(centered_X)
  }
  return(apply(centered_X,margin,function(x,N){
return((1/(sqrt((sum(x^2)/(N-1)))))*1/(sqrt(N-1)))
  },N=N))
}

alt_AF <- function(X){
  colSums(X)/(2*nrow(X))
  #return(apply(X,MARGIN = margin,))
}

alt_S <- function(X,margin=2){
  if(margin==2){
    N <- nrow(X)
  }else{
    N <- ncol(X)
  }
  return(apply(X,margin,function(x,N){
  tx <-sum(x)/(2*N)
  return(2*tx*(1-tx))
  },N=N))
}


gen_ty_h5 <- function(snp_df,snp_h5file,beta_h5file,tparam_df,AF=numeric(),sim_ind,chunksize=1000){

    p <- as.integer(nrow(snp_df))
    N <-as.integer(unique(tparam_df$n))
    g <- as.integer(nrow(tparam_df))
    stopifnot(min(sim_ind)>0,
              length(sim_ind)<=N,
              all(tparam_df$p==p),
              length(N)==1)
    snp_lff  <- BBmisc::chunk(snp_df$snp_id,chunk.size=chunksize) %>% purrr::map(~list(subset_rows=.x,
                                                                                filename=snp_h5file,
                                                                                datapath="dosage",
                                                                                subset_cols=sim_ind))
    beta_lff  <- BBmisc::chunk(1:p,chunk.size=chunksize) %>% purrr::map(~list(subset_cols=.x,
                                                                       filename=beta_h5file,
                                                                       datapath="Beta"))
    S_lff  <- BBmisc::chunk(1:p,chunk.size=chunksize) %>% purrr::map(~list(subset=.x,
                                                                       filename=beta_h5file,
                                                                       datapath="S"))

    stopifnot(all(purrr::map2_lgl(snp_lff,beta_lff,~length(.x$subset_rows)==length(.y$subset_cols))))
    stopifnot(all(purrr::map2_lgl(snp_lff,S_lff,~length(.x$subset_rows)==length(.y$subset))))

    dims_beta <-as.integer(c(g,p))

    EigenH5::create_matrix_h5(filename = beta_h5file,
                              groupname = "/",
                              dataname = "Beta",
                              data = numeric(),
                              dims = dims_beta,
                              chunksizes=as.integer(c(g,100)))
    EigenH5::write_vector_h5(filename = beta_h5file, groupname = "/", dataname="S",data = runif(p))
    ymat <- simulate_y_h5(list(SNP=snp_lff,Beta=beta_lff,S=S_lff),p,N,g,tparam_df$tsigu,Af=AF)
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


gen_sim_phenotype_h5 <- function(snp_df,snp_h5file,beta_h5file,tparam_df,AF=numeric(),sim_ind=(1:(tparam_df$n[1])),chunksize=1000){

    ty <- gen_ty_h5(snp_df,snp_h5file,beta_h5file,tparam_df,AF=AF,sim_ind=sim_ind,chunksize=chunksize)
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







