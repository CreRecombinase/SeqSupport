
#' Creates a generator for LDshrink
#' @param filename HDF5 filename
#' @param SNP_path datapath to haplotype data
#' @param snp_df dataframe specifying a subset of SNPs to compute LD for
LDshrink_h5_generator <- function(filename,SNP_subset = NULL,SNP_path="dosage", SNPinfo_path = "SNPinfo",map_file=filename,map_path="SNPinfo/map",mapinfo_path="SNPinfo", ...){
  subset_index <- gen_subset_index(SNP_subset,filename,SNPinfo_path,...)
}


#' Generate subset SNP
#' Generate an index for HDF5 to subset data
#' @param SNP_subset
#' @param filename path to HDF5 file we're going to subset
#' @param SNPinfo_path
#'
#' @return
#' @export
#'
#' @examples
gen_subset_index <-function(SNP_subset,filename, SNPinfo_path,SNP_path,...){
  UseMethod("gen_subset_index")
}

gen_subset_index.data.frame <- function(SNP_subset,filename,SNPinfo_path,...){
  snp_info_cols <- EigenH5::ls_h5(filename, SNPinfo_path)
  subset_cols <- colnames(SNP_subset)

  max_p <- EigenH5::get_dims_h5(filename,fs::path(SNPinfo_path,snp_info_cols[1]))
  if("snp_id" %in% subset_cols){
    stopifnot(all(dplyr::between(SNP_subset$snp_id,1,max_p)))
    return(SNP_subset$snp_id)
  }
  if(!("snp_id" %in% snp_info_cols)){
    snp_df <- dplyr::inner_join(SNP_subset, dplyr::mutate(EigenH5::read_df_h5(filename, SNPinfo_path),snp_id=1:n()))
  }else{
    snp_df <- dplyr::inner_join(SNP_subset, EigenH5::read_df_h5(filename, SNPinfo_path))
  }
  ncols <-colnames(snp_df)
  if(is.null(snp_df[["region_id"]])){
    snp_df <-  LDshrink::assign_snp_block(snp_df,break_df,assign_all=T)
  }
  if(is.null(snp_df[["map"]])){
    snp_df <- LDshrink::assign_map(snp_df,map_df)
  }
  dplyr::select(snp_df,region_id,selection_id=snp_id,map)
  # return(snp_df$snp_id)
}


ex_subset_function <- function(snp_df, AF_cutoff=0.05, ...){
  if(!is.null(snp_df[["selection_id"]])){
    dplyr::mutate(snp_df,AF=Allele_Freq) %>%
      dplyr::filter(AF>=AF_cutoff) %>%
      dplyr::select(region_id,selection_id,map) %>% return()
  }
  Allel_Freq <- snp_df[["AF"]]
  if(is.null(Allel_Freq)){
    dosage_mat <-  EigenH5::read_matrix_h5(filename,
                                           SNP_path,
                                           subset_rows=snp_df$snp_id)
    dos_dims <-dim(dosage_mat)
    if(SNPfirst){
      sum_func <-rowSums
      N <-dos_dims[2]
    }else{
      sum_func <- colSums
      N <- dos_dims[1]
    }
    Allele_Freq <- sum_func(dosage_mat)/(2*N)
  }
  if(is.null(snp_df[["region_id"]])){
    snp_df <-  LDshrink::assign_snp_block(snp_df,break_df,assign_all=T)
  }
  if(is.null(snp_df[["map"]])){
    snp_df <- LDshrink::assign_map(snp_df,map_df)
  }
  dplyr::mutate(snp_df,AF=Allele_Freq) %>%
    dplyr::filter(AF>=AF_cutoff) %>%
    dplyr::select(region_id,selection_id=snp_id,map)
}

gen_subset_index.function <- function(SNP_subset,filename,SNPinfo_path,SNP_path,...){
  snp_info_cols <- EigenH5::ls_h5(filename, SNPinfo_path)
  max_p <- EigenH5::get_dims_h5(filename,fs::path(SNPinfo_path,snp_info_cols[1]))
  SNPfirst <- EigenH5::get_dims_h5(filename,SNP_path)[1]==max_p
  chunksize <- 1000
  chunkl <- split(1:max_p,ggplot2::cut_interval(1:max_p,length=chunksize,labels=FALSE))
  pf <- purrr::partial(SNP_subset,filename = filename,SNPfirst = SNPfirst, SNP_path = SNP_path, .lazy = F,.first=F)
  res_l <- furrr::future_map(chunkl,
                             function(ind){
                               pf(snp_df = EigenH5::read_df_h5(filename, SNPinfo_path, subset=ind))
                             },
                             .options = furrr::future_options(packages="EigenH5"))
  return(unlist(res_l))
}

gen_subset_index.list <- function(SNP_subset,filename,SNPinfo_path,SNP_path,...){
  argl <- list(...)
  purrr::map(SNP_subset,
             ~purrr::invoke(gen_subset_index(filename=filename,SNPinfo_path=SNPinfo_path,SNP_path=SNP_path),SNP_subset=.x,
                            argl))
}

gen_subset_index.integer <-function(SNP_subset,filename,SNPinfo_path,SNP_path){
  snp_info_cols <- EigenH5::ls_h5(filename, SNPinfo_path)
  subset_cols <- colnames(SNP_subset)
  max_p <- EigenH5::get_dims_h5(filename,fs::path(SNPinfo_path,snp_info_cols[1]))
  stopifnot(dplyr::between(SNP_subset,1,max_p))
  snp_df <- read_df_h5(filename,SNPinfo_path,subset=SNP_subset)
  if(is.null(snp_df[["region_id"]])){
    snp_df <-  LDshrink::assign_snp_block(snp_df,break_df,assign_all=T)
  }
  if(is.null(snp_df[["map"]])){
    snp_df <- LDshrink::assign_map(snp_df,map_df)
  }
  dplyr::select(snp_df,region_id,selection_id=snp_id,map)
}

gen_subset_index.NULL <- function(SNP_subset,filename,SNPinfo_path, SNP_path){
  snp_info_cols <- EigenH5::ls_h5(filename, SNPinfo_path)
  subset_cols <- colnames(SNP_subset)
  max_p <- EigenH5::get_dims_h5(filename,fs::path(SNPinfo_path,snp_info_cols[1]))
  return(list(snp_id=1:max_p))
}

gen_subset_index.default <-function(SNP_subset,filename,SNPinfo_path, SNP_path){
  stop("Don't know how to generate subset_index for ",class(SNP_subset)[[1]],call. = FALSE)
}





offset_datasize_generator <- function(filename,generator_datapath,split_datapath,...){
  split_mat <- EigenH5::chunk_vec(EigenH5::read_vector_h5(filename,split_datapath))
  i <- 0
  num_chunks <- nrow(split_mat)
  function() {
    i <<- i + 1
    read_matrix_h5(filename,generator_datapath,offset_row=split_mat[i,1],datasize_row=split_mat)
  }
}
