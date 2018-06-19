#' read_SNPinfo
#' reads SNPinfo from an HDF5 file
#' @param snpfile
#' @param chr_to_char boolean specifying whether to convert chrom to character type (by prefixing "chr"),default is TRUE
#' @param extra_cols character vector specifying other columns you'd like returned
#' @return dataframe with (at least) `chr`,`pos` and optionally `extra_cols`
read_SNPinfo <- function(snpfile,chr_to_char=T, extra_cols = NULL, id_col=NULL){
  library(dplyr)
  # snpfile_dsets <- RcppEigenH5::h5ls(snpfile)

  pos <- read_ivec(snpfile,"/","pos")
  chr <- read_ivec(snpfile,"/","chr")
  snp_df <- tibble::data_frame(pos = pos, chr = chr) %>% dplyr::mutate(snp_id = paste0("SNP: ",1:n()))
  if (chr_to_char) {
    snp_df <- dplyr::mutate(snp_df,chr = paste0("chr",chr))
  }
  return(snp_df)
}




# mach2h5 <- function(machfile.gz,hdf5file,datapath,subset_rows=NULL,subset_cols=NULL){
#
#
# read_delim_chunked(machfile.gz,)
#
# }


waitr <- function(listres,timeout=Inf){
  num_p <- length(listres)
  num_f <-sum(purrr::map_lgl(listres,resolved))
  pb <- progress::progress_bar$new(total = num_p)
  while(num_f<num_p){
    Sys.sleep(0.5)
    num_f <-sum(purrr::map_lgl(listres,resolved))
    pb$update(num_f/num_p)
  }
  return(listres)
}


read_SNPinfo_allel <- function(snpfile,groupname="variants"){

  objd <- EigenH5::get_objs_h5(snpfile,"variants")
  obj_dims <- purrr::map_int(objd,~length(EigenH5::get_dims_h5(snpfile,"variants",.x)))
  obj_t <- purrr::map_int(objd,~as.integer(EigenH5::check_dtype(snpfile,"variants",.x)))
  ALT <- read_vector_h5(snpfile,"variants","ALT")
  all_t <- objd[obj_dims==1&(obj_t %in%c(13,14,16))]
  known_weird_cols <- c("AA","AN","CS","DP","END","MC","MEND","MLEN","MSTART","NS","QUAL","SVLEN","SVTYPE","TSD")
  snp_df <-EigenH5::read_df_h5(snpfile,"variants",subcols=all_t)%>%
    dplyr::rename(SNP=ID,chr=CHROM,pos=POS) %>%
    dplyr::select(-dplyr::one_of(known_weird_cols)) %>%
    tidyr::unite(allele,REF,ALT,sep=",") %>% mutate(snp_id=1:n(),chr=as.integer(chr))

  return(snp_df)
}



match_df_h5 <- function(df,output_filenames,
                        output_group_prefix,
                        output_datanames,isVector=F){
    chunksize_max<-1024
    output_dff <-dplyr::mutate(df,
                               create_dynamic=F,
                               row_offsets=0L,
                               col_offsets=0L,
                               datatypes="numeric",
                               col_chunksizes=row_chunksizes,
                               col_c_chunksizes=pmin(chunksize_max,col_chunksizes),
                               row_c_chunksizes=col_c_chunksizes,
                               filenames=output_filenames,
                               groupnames=paste0(output_group_prefix,"/",chunk_group),
                               datanames=output_datanames)
    if(isVector){
        output_dff <- output_dff %>% dplyr::mutate(col_chunksizes=1L,col_c_chunksizes=1L)
    }
    return(output_dff)
}



cross_df_h5 <- function(df_rows,df_cols,output_filenames,
                     output_groupnames,
                     output_datanames){
    uh_dff <- dplyr::inner_join(dplyr::mutate(df_cols,col=NA,mid=1:n()),
                                dplyr::mutate(df_rows,col=NA,mid=1:n()),
                                by="col",suffix=c("_cols","_rows")) %>%
        dplyr::mutate(filenames=output_filenames,
                      groupnames=output_groupnames,
                      datanames=output_datanames)

    ret_dff <- uh_dff %>% dplyr::arrange(mid_cols,mid_rows)%>%
      dplyr::select(filenames,groupnames,datanames,
                                 row_offsets=row_offsets_rows,
                                 row_chunksizes=row_chunksizes_rows,
                                 col_offsets=row_offsets_cols,
                                 col_chunksizes=row_chunksizes_cols)
    return(ret_dff)
}


read_vec <- function(h5filename,datapath){
  if(substr(datapath,1,1)!="/"){
    datapath <- paste0("/",datapath)
  }
  return(c(rhdf5::h5read(h5filename,datapath)))
}


prep_h5file <- function(h5filename,create_dir=F){

    if(class(h5filename)[1]=="character"){


        if(!file.exists(h5filename)){
            file_dir <- dirname(h5filename)
            if(!dir.exists(file_dir)){
                if(!create_dir){
                    stop(paste0("Directory ",file_dir," does not exist, set create_dir=T"))
                }else{
                    dir.create(file_dir,recursive = T)
                }
            }
            rhdf5::h5createFile(h5filename)
        }else{
            if(!rhdf5::H5Fis_hdf5(h5filename)){
                file.remove(h5filename)
                rhdf5::h5createFile(h5filename)
            }
        }
        return(T)
    }
    if(class(h5filename)[1]=="H5IDComponent"){
        return(T)
    }
}

get_h5f <- function(hdf5file,readOnly=T){
  read_wr <- ifelse(readOnly,"H5F_ACC_RDONLY","H5F_ACC_RDWR")
  if(class(hdf5file)[1]=="character"){
    stopifnot(file.exists(hdf5file))
    return(rhdf5::H5Fopen(hdf5file,read_wr))
  }else{
    stopifnot(class(hdf5file)[1]=="H5IdComponent")
    return(hdf5file)
  }
}


is_transpose_h5 <- function(hdf5file,datapath){
  if(class(hdf5file)[1]=="character"){
    hf <- rhdf5::H5Fopen(hdf5file,"H5F_ACC_RDONLY")
  }else{
    hf <- hdf5file
  }
  if(class(datapath)[1]=="character"){
    hd <- rhdf5::H5Dopen(hf,datapath)
  }else{
    hd <- datapath
  }
  attr <- rhdf5::H5Aopen(hd,"doTranspose")
  doTranspose <- rhdf5::H5Aread(attr)==0
  rhdf5::H5Aclose(attr)
  if(class(datapath)[1]=="character"){
    rhdf5::H5Dclose(hd)
  }
  if(class(hdf5file)[1]=="character"){
    rhdf5::H5Fclose(hf)
  }
  return(doTranspose)
}

get_rownum_h5 <- function(hdf5file,datapath){
  if(class(hdf5file)[1]=="character"){
    hf <- rhdf5::H5Fopen(hdf5file,"H5F_ACC_RDONLY")
  }else{
    hf <- hdf5file
  }
  if(class(datapath)[1]=="character"){
    hd <- rhdf5::H5Dopen(hf,datapath)
  }else{
    hd <- datapath
  }

  hd <- rhdf5::H5Dopen(hf,datapath)
  hs <- rhdf5::H5Dget_space(hd)
  dim <- rhdf5::H5Sget_simple_extent_dims(hs)$size
  if(is_transpose(hf,hd)){
    retdim <- dim[2]
  }else{
    retdim <- dim[1]
  }
  if(class(datapath)[1]=="character"){
    rhdf5::H5Dclose(hd)
  }
  if(class(hdf5file)[1]=="character"){
    rhdf5::H5Fclose(hf)
  }
  return(retdim)

}

get_colnum_h5 <- function(hdf5file,datapath){
  if(class(hdf5file)[1]=="character"){
    hf <- rhdf5::H5Fopen(hdf5file,"H5F_ACC_RDONLY")
  }else{
    hf <- hdf5file
  }
  if(class(datapath)[1]=="character"){
    hd <- rhdf5::H5Dopen(hf,datapath)
  }else{
    hd <- datapath
  }

  hd <- rhdf5::H5Dopen(hf,datapath)
  hs <- rhdf5::H5Dget_space(hd)
  dim <- rhdf5::H5Sget_simple_extent_dims(hs)$size
  if(is_transpose(hf,hd)){
    retdim <- dim[1]
  }else{
    retdim <- dim[2]
  }
  if(class(datapath)[1]=="character"){
    rhdf5::H5Dclose(hd)
  }
  if(class(hdf5file)[1]=="character"){
    rhdf5::H5Fclose(hf)
  }
  return(retdim)

}






split_i <- function(p,i=1,chunksize,retl=list()){
  if(i+chunksize-1>=p){
    retl[[length(retl)+1]] <- i:p
    return(retl)
  }
  retl[[length(retl)+1]] <- i:(i+chunksize-1)
  return(split_i(p,i+chunksize,chunksize,retl))
}



#
# write_vec <- function(h5filename,groupname="/",dataname,data,deflate_level=0L,create_dir=F){
#   library(rhdf5)
#     prep_h5file(h5filename,create_dir)
#
#     groupname <- ifelse(groupname=="/","",groupname)
#     d_path <- paste0(groupname,"/",dataname)
#     if(!group_exists(h5filename,groupname)){
#       rhdf5::h5createGroup(h5filename,groupname)
#     }
#   groupname <- ifelse(groupname=="/","",groupname)
#   d_path <- paste0(groupname,"/",dataname)
#   if(typeof(data)=="character"){
#     rhdf5::h5createDataset(
#       file=h5filename,
#       dataset=d_path,showWarnings=F,
#       dims=as.integer(length(data)),maxdims=as.integer(length(data)),
#       storage.mode=storage.mode(data),chunk = as.integer(length(data)/10),size=255L,level=deflate_level)
#   }else{
#     rhdf5::h5createDataset(
#       file=h5filename,
#       dataset=d_path,showWarnings=F,
#       dims=as.integer(length(data)),maxdims=as.integer(length(data)),
#       storage.mode=storage.mode(data),chunk = as.integer(length(data)/10),level=deflate_level)
#   }
#   rhdf5::h5write(data,file=h5filename,name=d_path)
# }




group_exists <- function(h5file,groupname){
  if(groupname=="/"){
    return(T)
  }

  if(substr(groupname,1,1)!="/"){
    groupname <- paste0("/",groupname)
  }
  h5gs <- dplyr::filter(rhdf5::h5ls(h5file),otype=="H5I_GROUP") %>%
    dplyr::mutate(group=ifelse(group=="/",group,paste0(group,"/")))
  # if(sum(grepl(groupname,))

  h5g <- paste0(h5gs$group,h5gs$name)
  return(any(h5g==groupname))
}





chunk_df_h5 <- function(filename,groupname,dataname,chunksize_row=NULL,chunksize_col=NULL){
  data_dims <- EigenH5::get_dims_h5(filename = filename,groupname = groupname,dataname = dataname)
  chunksize_row <- min(c(chunksize_row,data_dims[1]))
  is_mat <- length(data_dims)==2
  if(is_mat){
    chunksize_col <- min(c(chunksize_col,data_dims[2]))
  }
  row_offset <-as.integer(seq.int(from = 0,to = data_dims[1]-1,by=chunksize_row))
  row_chunksize <-as.integer(pmin(chunksize_row,data_dims[1]-row_offset))
  h_df <-tibble::data_frame(row_offsets=row_offset,
                              row_chunksizes=row_chunksize,
                              filenames=filename,
                              groupnames=groupname,
                              datanames=dataname)

  if(is_mat){
    col_offset <-as.integer(seq.int(from = 0,to = data_dims[2]-1,by=chunksize_col))
    col_chunksize <-as.integer(pmin(chunksize_col,data_dims[2]-col_offset))
    h_df <-tibble::data_frame(col_offsets=col_offset,
                                col_chunksizes=col_chunksize,
                                filenames=filename,
                                groupnames=groupname,
                                datanames=dataname) %>%
      dplyr::inner_join(h_df)
  }
  return(h_df)
}


gen_outer_df <-function(col_chunk_df,row_chunk_df){
  dplyr::inner_join(dplyr::select(col_chunk_df,col_offsets,col_chunksizes) %>% dplyr::mutate(c=NA),dplyr::select(row_chunk_df,row_offsets,row_chunksizes)%>% dplyr::mutate(c=NA)) %>%
  dplyr::select(-c) %>% return()
}


gen_map_eqtl_df <- function(snpfile,expfile,uhfile,snp_chunksize=50000,exp_chunksize=10000){

  snp_dims <-EigenH5::get_dims_h5(filename = snpfile,"/","dosage")
  exp_dims <-EigenH5::get_dims_h5(expfile,"trait","ymat")
  stopifnot(exp_dims[1]==snp_dims[2])
  p <-snp_dims[1]
  g <- exp_dims[2]
  snp_chunksize <- min(snp_chunksize,p)
  exp_chunksize <- min(exp_chunksize,g)



}


#
# read_df_h5 <- function(h5filepath,groupname=NULL,subcols=NULL,filtervec=list(NULL)){
#   library(rhdf5)
#   stopifnot(file.exists(h5filepath))
#   if(is.null(groupname)){
#     groupname <- "/"
#   }
#   if(substr(groupname,1,1)!="/"){
#     groupname <- paste0("/",groupname)
#   }
#   stopifnot(group_exists(h5filepath,groupname))
#
#   dsets <- list.datasets(h5filepath,groupname,subcols)
#   if(substr(groupname,nchar(groupname),nchar(groupname))!="/"){
#     groupname <- paste0(groupname,"/")
#   }
#   paths <- as.list(paste0(groupname,dsets))
#   names(paths) <- dsets
#   stopifnot(length(dsets)>0)
#   if(is.logical(filtervec)){
#     filtervec <- purrr::map(filtervec,which)
#   }
#   return(purrr::map_dfc(.x = paths,purrr::compose(c,rhdf5::h5read),file=h5filepath,index=filtervec))
# }


## #' write_list_h5
## #' write a list to an HDF5 file
## #' @param h5file the name of the HDF5 file to write
## #' @param groupname the name of the group for the HDF5 file (can specify several groups, e.g 'a/b/c')
## #' @param dataname the name of the matrix in the file
## #' @param data the matrix to write, can be of any atomic type
## #' @param deflate_level integer specifying compression level
## #' @param doTranspose Bool indicating whether or not to write the data in row-major or column major order
## write_list_h5 <- function(h5file, groupname="/", dataname, data, deflate_level = as.integer(c(0))){





## }



#'
#' #' write_mat_h5
#' #' write a matrix to an HDF5 file
#' #' @param h5file the name of the HDF5 file to write
#' #' @param groupname the name of the group for the HDF5 file (can specify several groups, e.g 'a/b/c')
#' #' @param dataname the name of the matrix in the file
#' #' @param data the matrix to write, can be of any atomic type
#' #' @param deflate_level integer specifying compression level
#' #' @param doTranspose Bool indicating whether or not to write the data in row-major or column major order
#' write_mat_h5 <- function(h5file, groupname="/", dataname, data, deflate_level = as.integer(c(0)),doTranspose=T){
#'     library(rhdf5)
#'     prep_h5file(h5file,create_dir=T)
#'
#'
#'     groupname <- ifelse(groupname=="/","",groupname)
#'     if(!group_exists(h5file,groupname)){
#'         rhdf5::h5createGroup(h5file,groupname)
#'     }
#'     d_path <- paste0(groupname,"/",dataname)
#'     if(doTranspose){
#'         rhdf5::h5write(data,h5file,d_path)
#'     }else{
#'         rhdf5::h5write(t(data),h5file,d_path)
#'     }
#'
#'     fh <- rhdf5::H5Fopen(h5file,flags = "H5F_ACC_RDWR")
#'     dsp <-rhdf5::H5Oopen(fh,d_path)
#'     transpose_wr <- ifelse(doTranspose,1L,0L)
#'     rhdf5::h5writeAttribute(transpose_wr,dsp,"doTranspose")
#'     rhdf5::H5Oclose(dsp)
#'     rhdf5::H5Fclose(fh)
#' }
#'
#' #' write_df_h5
#' #' Write a dataframe to an HDF5 file
#' #' @param df dataframe to write (current support for nested/list dataframes is iffy)
#' #' @param groupname name of group in which to write the dataframe, this will basically be the name of the dataframe in the HDF5 file
#' #' @param outfile the name of the hdf5 file to write
#' #' @param deflate_level integer specifying level of compression to apply, with higher indicating higher compression
#' write_df_h5 <- function(df,groupname="/",outfile,deflate_level=4L){
#'
#'   # prep_h5file(outfile,create_dir = T)
#'
#'
#'   # dataname <- colnames(df)
#'   # datapaths <-paste0(groupname,"/",dataname)
#'
#'   purrr::iwalk(df,
#'                function(val,name,filename,groupname,deflate_level){
#'                  EigenH5::write_vector_h5(filename = filename,
#'                                           groupname = groupname,
#'                                           dataname = name,
#'                                           data = val)
#'                },
#'                filename=outfile,
#'                groupname=groupname,
#'                deflate_level=deflate_level)
#' }


gds2hdf5 <- function(gdsfile,hdf5file,deflate_level=4L){
  # library(rhdf5)
  is_c <- class(gdsfile)=="character"
  if(is_c){
    gds <- SeqArray::seqOpen(gdsfile)
  }else{
    gds <- gdsfile
  }
    snp_info <-read_SNPinfo_gds(gds,alleles=T,MAF=T,region_id=F,map = F,info=T,more=list(rs="annotation/id")) %>%
        dplyr::mutate(chr=as.integer(chr)) %>% dplyr::mutate(snp_id=as.integer(ifelse(!is.integer(snp_id),1:n(),snp_id))) %>%
        dplyr::arrange(chr,pos) %>% dplyr::mutate(nsnp_id=1:n())
    stopifnot(dplyr::group_by(snp_info,chr) %>%
              dplyr::summarise(is_sorted=!is.unsorted(snp_id)) %>%
              dplyr::summarise(is_sorted=all(is_sorted)) %>%
              dplyr::pull(1))
    dosage2hdf5(gds=gds,hdf5file=hdf5file,snp_info=snp_info)
    snp_info <- dplyr::mutate(snp_info,snp_id=nsnp_id) %>% dplyr::select(-nsnp_id)
    sample_info <- tibble::data_frame(sample_id=seqGetData(gds,"sample.id"))
    EigenH5::write_df_h5(sample_info,"SampleInfo",hdf5file)
    EigenH5::write_df_h5(df = snp_info,groupname = "SNPinfo",hdf5file)
    #Close it if we opened it
    if(is_c){
      seqClose(gds)
    }
}


read_quh_df <- function(h5file,subvecs=NULL,pvv=1){

  tparam_df <- EigenH5::read_df_h5(h5file,"SimulationInfo") %>% group_by(tpve) %>% mutate(pve_replicate=1:n()) %>% ungroup()
  if(is.null(subvecs)){
    subvecs <- 1:nrow(tparam_df)
  }
  exec_fn <-  function(quhl,tpdfl,D,region_id){
    data_frame(fgeneid=tpdfl$fgeneid,
               region_id=region_id,
               quh=quhl,
               D=D,
               p_n=tpdfl$p/tpdfl$n,
               tbias=tpdfl$tbias,
               pve_replicate=tpdfl$pve_replicate,
               tpve=tpdfl$tpve,
               tsigu=tpdfl$tsigu)%>%dplyr::mutate(eigen_id=1:n())%>%
      RSSp::filter_pvv(pvv) %>%return()
  }
  purrr::map2_df(
    purrr::array_branch(EigenH5::read_matrix_h5(h5file,"/","quh",subset_cols = subvecs),margin=2),
    purrr::transpose(dplyr::slice(tparam_df,subvecs)),
    exec_fn,
    D=read_vector_h5(h5file,"/","D"),
    region_id=read_vector_h5(h5file,"SNPinfo","region_id")) %>% return()

}



read_uh_df <- function(h5file,subvecs=NULL,evdf=NULL){

  tparam_df <- EigenH5::read_df_h5(h5file,"SimulationInfo") %>% group_by(tpve) %>% mutate(pve_replicate=1:n()) %>% ungroup()
  if(is.null(subvecs)){
    subvecs <- 1:nrow(tparam_df)
  }
  exec_fn <-  function(uhl,tpdfl,L2=NULL,region_id){
    data_frame(fgeneid=tpdfl$fgeneid,
               uh=uhl,
               L2=L2,
               p_n=tpdfl$p/tpdfl$n,
               pve_replicate=tpdfl$pve_replicate,
               tbias=tpdfl$tbias,
               tpve=tpdfl$tpve,
               tsigu=tpdfl$tsigu) %>% dplyr::mutate(snp_id=1:n())%>% return()
  }
    L2 <- NULL
    if(!is.null(evdf)){
        L2 <- purrr::as_vector(
                         purrr::map(EigenH5::ls_h5(evdf,"L2",full_names=T),~EigenH5::read_vector_h5(evdf,.x,"L2"))
                     )
        }
  purrr::map2_df(
    purrr::array_branch(EigenH5::read_matrix_h5(h5file,"/","uh",subset_cols = subvecs),margin=2),
    purrr::transpose(dplyr::slice(tparam_df,subvecs)),
    exec_fn,
    L2=L2) %>% return()
}

dosage2hdf5 <- function(gds,hdf5file,chunksize=c(150),snp_info){
  #library(rhdf5)
  p <- calc_p(gds)
  N <- calc_N(gds)
  dims <- as.integer(c(p,N))
  is_haplo <- is_haplo(gds)
  if(!dir.exists(dirname(hdf5file))){
    dir.create(dirname(hdf5file))
  }
  if(!file.exists(hdf5file)){
   # h5createFile(hdf5file)
  }
  chunksize_max <- 1024^2
  if(prod(c(chunksize,N))>chunksize_max){
    tchunksize <-as.integer(ceiling((chunksize_max/N)/2))
    warning(paste0("chunksize: ",chunksize," is too big, resizing to",tchunksize))
    chunksize <-tchunksize
  }
  EigenH5::create_matrix_h5(filename = hdf5file,
                            groupname = "/",
                            dataname = "dosage",
                            dims = dims,
                            data=numeric(),
                            chunksizes = as.integer(c(chunksize,N)))
  write_chunk <-function(index,x,h5loc,is_haplo,N){
    if(is_haplo){
      tobj <- t(2.0-x)
    }else{
      tobj <- t(x)
    }
    tp <-nrow(tobj)
    oindex_r <-snp_info$nsnp_id[seq.int(from=index,length.out = tp)]
    oindex <- oindex_r[1]
    stopifnot(all(oindex_r==seq.int(from=oindex,length.out = tp)))
    stopifnot(ncol(tobj)==N)
    # stopifnot(all(tobj>=0))
    EigenH5::write_matrix_h5(filename = h5loc,
                                  groupname="/",
                                  dataname="dosage",
                                  offsets = c(oindex-1,0L),
                             data = tobj)

  }
  # chrs <-as.character(unique(snp_info$chr))

  seqBlockApply(gdsfile = gds,
                var.name = "$dosage",
                FUN =write_chunk,
                h5loc=hdf5file,
                is_haplo=is_haplo,
                N=N,
                var.index="relative")

}







