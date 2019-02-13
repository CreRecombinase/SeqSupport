map_eqtl_h5 <- function(snp_h5,exp_h5,
                        covar_h5,uh_h5,
                        ncovar=0,
                        snp_path="dosage",
                        exp_path="trait/ymat",
                        covar_path="covariates",
                        snp_info="SNPinfo",
                        exp_info="TraitInfo",
                        snp_chunksize=5000,
                        exp_chunksize=Inf,
                        subgwasf=NULL,
                        subsnpf=NULL,
                        progress=TRUE,
                        threads=parallel::detectCores(),
                        ...){



    exp_df <-read_df_h5(exp_h5,exp_info)
    write_df_h5(exp_df, uh_h5, "Traitinfo")
    g <- nrow(exp_df)

    expdims <- dim_h5(exp_h5,exp_path)

    EXPfirst <- g==expdims[1]
    if(EXPfirst){
        N  <- expdims[2]
    }else{
        N <- expdims[1]
    }


    if(ncovar>0){
        cov_dim <- dim_h5(covar_h5,covar_path)
        covarFirst <- cov_dim[2]==N
        if(covarFirst){
            covar_tot <- cov_dim[1]
            stopifnot(ncovar<=covar_tot)
            covarmat <- t(read_matrix_h5(covarf,covar_path,subset_rows=seq_len(ncovar)))
        }else{
            covar_tot <- cov_dim[2]
            stopifnot(ncovar<=covar_tot)
            covarmat <- read_matrix_h5(covarf,covar_path,subset_cols=seq_len(ncovar))
        }
    }else{
        covarmat <- matrix(1,N,1)
    }

    if(!is.null(subgwasf)){
        if(tools::file_ext(subgwasf)=="h5"){
            ind_v <- read_vector_h5(subgwasf,"SampleInfo/sample_id")
        }else{
            ind_v <- readRDS(subgwasf)
        }
    }else{
        ind_v  <- 1:N
    }

    if(!is.null(subsnpf)){
        if(tools::file_ext(subsnpf)=="h5"){
            snp_df <- read_df_h5(subsnpf,snp_info)
        }else{
            snp_df <- read_delim(subsnpf,delim="\t")
        }
    }else{
        snp_df <- EigenH5::read_df_h5(snp_h5,snp_info)
    }

    p <- nrow(snp_df)
    dim_dosage <- EigenH5::dim_h5(snp_h5,snp_path)
    SNPfirst <-dim_dosage[2]==N


    snp_chunks <- max(c(ceiling(p/snp_chunksize),1))
    exp_chunks  <-max(c(ceiling(g/exp_chunksize),1))
    write_df_h5(snp_df,uh_h5,"SNPinfo")


    snp_df <- snp_df %>%  dplyr::mutate(snp_chunk=ggplot2::cut_number(x=snp_id,n=snp_chunks,labels=F),snp_chunk_id=1:n())
    snp_l <-  split(select(snp_df,snp_id,snp_chunk_id),snp_df$snp_chunk)

    if(SNPfirst){
        snp_lff  <- snp_l %>% purrr::map(~list(subset_rows=.x$snp_id,
                                               filename=snp_h5,
                                               datapath=snp_path))
    }else{
        snp_lff  <- snp_l %>% purrr::map(~list(subset_cols=.x$snp_id,
                                               filename=snp_h5,
                                               datapath=snp_path))
    }
    if(nrow(exp_df)>1){
      exp_df <- exp_df %>%  dplyr::mutate(exp_chunk=ggplot2::cut_number(trait_id,exp_chunks,labels=F),exp_chunk_id=1:n())
    }else{
      exp_df <- dplyr::mutate(exp_df,exp_chunk=1L,exp_chunk_id=1L)
    }
    exp_l <-  split(select(exp_df,trait_id,exp_chunk,exp_chunk_id),exp_df$exp_chunk)

    if(EXPfirst){
      exp_lff  <- exp_l %>% purrr::map(~list(subset_rows=.x$trait_id,
                                             filename=exp_h5,
                                             datapath=exp_path))
    }else{
      exp_lff  <- exp_l %>% purrr::map(~list(subset_cols=.x$trait_id,
                                             filename=exp_h5,
                                             datapath=exp_path))
    }


    snp_exp_lff <-purrr::cross2(snp_l,exp_l)


    uh_lff <- map(snp_exp_lff,~list(subset_rows=.x[[1]]$snp_chunk_id,subset_cols=.x[[2]]$exp_chunk_id,filename=uh_h5,datapath="uh"))
    se_lff <-map(uh_lff,~update_list(.,datapath="se"))
    EigenH5::create_matrix_h5(uh_h5,"uh", numeric(),dims=as.integer(c(p,g)))
    EigenH5::create_matrix_h5(uh_h5,"se", numeric(),dims=as.integer(c(p,g)))
    map_eQTL_chunk_h5(snp_lff,exp_lff,uh_lff,se_lff,covarmat,
                      options=list(SNPfirst=SNPfirst,
                                   EXPfirst=EXPfirst,
                                   progress=progress,
                                   threads=threads))


}
