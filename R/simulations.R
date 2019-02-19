#' Generate parameters of simulation
#' @param pve true PVE
#' @param bias true value for bias/confounding (as a percentage of true variance in normalized effect size)
#' @param nreps number of repetitions of each combo of parameters
#' @param n sample size
#' @param p number of SNPs
#' @param fgeneid name for each trait (defaults to 1:ntraits)
tparamdf_norm <- function(pve=NULL, bias = 0, nreps, n, p,pcovar=0,fgeneid=NULL,sigu=NULL){

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
        dplyr::mutate(fgeneid=as.character(1:n()), tsigu=sqrt(n/p*tpve)) %>%
        dplyr::select(-replicate) %>% mutate(pcovar=pcovar)
    if(!is.null(tfgeneid)){
        stopifnot(nrow(rfp)==length(tfgeneid))
        rfp <- dplyr::mutate(rfp,fgeneid=as.character(tfgeneid))
    }
    return(rfp)
}


#' calculate the scale of residuals, given X*Beta and the target PVE
#' @param vy vector with the variance of y (length is equal to the number of traits)
#' @param PVE vector specifying the target PVE (length should be the same as the number of traits)
gen_ti <- function(vy, PVE,vconf=rep(0L,length(vy))){
    stopifnot(length(vy)==length(PVE))
    # stopifnot(all(vconf[PVE==0]==0))
    retvec <- vy*(1/PVE-1)-vconf
    stopifnot(all(retvec[PVE>0]>=0))
    # retvec <- /vy-vy-vconf
    retvec[PVE==0] <- 1-vconf
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


gen_ty_h5 <- function(snp_df,snp_h5file,beta_h5file,tparam_df,sim_ind,chunksize=1000){

    p <- as.integer(nrow(snp_df))
    N <-as.integer(unique(tparam_df$n))
    g <- as.integer(nrow(tparam_df))
    num_chunks <- ceiling(p/chunksize)
    if(is.null(snp_df[["region_id"]])){
        snp_df <- mutate(snp_df, region_id = as.integer(gl(n=num_chunks,k=chunksize,length=p)))
    }
    stopifnot(exprs = {
      min(sim_ind) > 0
      length(sim_ind) <= N
      length(N) == 1
    })


    snp_l  <- split(snp_df$snp_id, snp_df$region_id) #%>% purrr::map(~list(subset_rows=.x,
    #                                                                                filename=snp_h5file,
    #                                                                               datapath="dosage",
    #                                                                              subset_cols=sim_ind))
    beta_l  <- split(1:p, snp_df$region_id) #%>% purrr::map(~list(subset_cols=.x,
    #         filename=beta_h5file,
    #                datapath="Beta"))

    betac_l  <- split(1:p, snp_df$region_id) #%>% purrr::map(~list(subset_cols=.x,
    #                         filename=beta_h5file,
    #                         datapath="Betac"))

    # D_l  <- paste0("EVD/",unique(snp_df$region_id),"/D")
    # #                       datapath=paste0("EVD/",.x,"/D")))
    #
    # Q_l  <- paste0("EVD/",unique(snp_df$region_id),"/Q")
    #%>% purrr::map_chr(~list(filename=evd_h5file,
    #                       datapath=paste0("EVD/",.x,"/Q")))


    S_l  <- split(1:p, snp_df$region_id)
    #%>% purrr::map(~list(subset=.x,
    #                      filename=beta_h5file,
    #                     datapath="S"))

    #stopifnot(all(purrr::map2_lgl(snp_lff,beta_lff,~length(.x$subset_rows)==length(.y$subset_cols))))
    #stopifnot(all(purrr::map2_lgl(snp_lff,S_lff,~length(.x$subset_rows)==length(.y$subset))))

    dims_beta <-as.integer(c(g,p))


    EigenH5::create_matrix_h5(filename = beta_h5file,
                              datapath="Beta",
                              data = numeric(),
                              dims = dims_beta,
                              chunksizes=as.integer(c(g,min(p,100))))
    # EigenH5::create_matrix_h5(filename = beta_h5file,
    #                           datapath="Betac",
    #                           data = numeric(),
    #                           dims = dims_beta,
    #                           chunksizes=as.integer(c(g,min(p,100))))

    EigenH5::write_vector_h5(data = runif(p),filename = beta_h5file, datapath="S")
    # ymat <- simulate_y_h5(list(SNP=snp_lff,
    #                            Beta=beta_lff,
    #                            Betac=betac_lff,
    #                            S=S_lff,
    #                            D=D_lff,
    #                            Q=Q_lff),p,N,g,tparam_df$tsigu,tsigc=tparam_df$tbias)
    cur_p_offset <- 0
    ymat <- matrix(0,n,g)
    pb <- progress_bar$new(total=length(snp_l))
    for(i in 1:length(snp_l)){
      X <- read_matrix_h5(snp_h5file,"dosage",subset_rows=snp_l[[i]],doTranspose=T)
      S <- (1/sqrt(colMeans(X^2)))*(1/sqrt(n-1))
      write_vector_h5(S,filename = beta_h5file,datapath = "S",subset=S_l[[i]])

      tp <- length(S)
      U <- matrix(rnorm(n = tp*g,mean=0,sd=rep(tparam_df$tsigu,times=p)),nrow = tp,ncol = g)
      Beta <- U*S
      # Betac <-matrix(rnorm(n = tp*g,mean=0,sd=rep(tparam_df$tbias,times=p)),nrow = tp,ncol = g)
      # Q <- read_matrix_h5(evd_h5file,datapath = Q_l[i])
      # D <- read_vector_h5(evd_h5file,datapath=D_l[i])
      # ymat <- ymat+X%*%(Beta+Q%*%(D*Betac))
      ymat <- ymat+X%*%Beta
      write_matrix_h5(t(Beta),filename = beta_h5file,"Beta",subset_cols=beta_l[[i]])
      # write_matrix_h5(t(Betac),filename = beta_h5file,"Betac",subset_cols=betac_l[[i]])
      pb$tick()
    }
    return(ymat)
}





gen_sim_resid_covar <- function(ty,tparam_df,C=matrix(1,nrow=nrow(ty),ncol=1),covar_wf = NULL){
  ncovar <- ncol(C)
  stopifnot(ncovar==1L)
  vy <- apply(ty, 2, var)
  n <- nrow(ty)
  stopifnot(ncol(ty)==nrow(tparam_df))
  residvec <- gen_ti(vy, tparam_df$tpve)
  if(is.null(tparam_df$covarp)){
    tparam_df$covarp <- 0
  }
  covarvec <- residvec*(tparam_df$covarp)
  covar_effect <-sapply(covarvec,function(ti){C*ti})
  if(!is.null(covar_wf)){
    covar_wf(covar_effect)
  }
  residvec <-(1-tparam_df$covarp)*residvec
  residmat <- sapply(residvec, function(ti, n){rnorm(n=n, mean=0, sd=sqrt(ti))}, n=n)
  ymat <- scale(ty+residmat, center=T, scale=F)
  return(ymat)
}


gen_sim_resid <- function(ty,tparam_df,Q=matrix(0,nrow=nrow(ty),ncol=1)){
    vy <- apply(ty, 2, var)
    N <- nrow(ty)
    stopifnot(ncol(ty)==nrow(tparam_df))
    K <- ncol(Q)
    Q_beta <- t(matrix(sqrt(((tparam_df$tbias)*(N-1))/K),
                       nrow=length(tparam_df$tbias),
                       ncol=K))

    Q_y <- Q%*%Q_beta


    residvec <- gen_ti(vy, tparam_df$tpve,vconf=tparam_df$tbias)



    residmat <- sapply(residvec, function(ti, n){rnorm(n=n, mean=0, sd=sqrt(ti))}, n=n)
    ymat <- scale(ty+residmat+Q_y, center=T, scale=F)
    return(ymat)
}


sim_phenotype_h5 <- function(snp_df,snp_h5file,beta_h5file,tparam_df,sim_ind=(1:(tparam_df$n[1])),chunksize=1000){

    ty <- gen_ty_h5(snp_df,snp_h5file,beta_h5file,tparam_df,sim_ind=sim_ind,chunksize=chunksize)

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







