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
  snp_df <- data_frame(pos = pos, chr = chr) %>% mutate(snp_id = paste0("SNP: ",1:n()))
  if (chr_to_char) {
    snp_df <- mutate(snp_df,chr = paste0("chr",chr))
  }
  return(snp_df)
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




write_vec <- function(h5filename,groupname="/",dataname,data,deflate_level=0L,create_dir=F){
  library(rhdf5)
    prep_h5file(h5filename,create_dir)

    groupname <- ifelse(groupname=="/","",groupname)
    d_path <- paste0(groupname,"/",dataname)
    if(!group_exists(h5filename,groupname)){
      rhdf5::h5createGroup(h5filename,groupname)
    }
  groupname <- ifelse(groupname=="/","",groupname)
  d_path <- paste0(groupname,"/",dataname)
  if(typeof(data)=="character"){
    rhdf5::h5createDataset(
      file=h5filename,
      dataset=d_path,showWarnings=F,
      dims=as.integer(length(data)),maxdims=as.integer(length(data)),
      storage.mode=storage.mode(data),chunk = as.integer(length(data)/10),size=255L,level=deflate_level)
  }else{
    rhdf5::h5createDataset(
      file=h5filename,
      dataset=d_path,showWarnings=F,
      dims=as.integer(length(data)),maxdims=as.integer(length(data)),
      storage.mode=storage.mode(data),chunk = as.integer(length(data)/10),level=deflate_level)
  }
  rhdf5::h5write(data,file=h5filename,name=d_path)
}

read_2d_mat_h5 <- function(h5filename,groupname="/",dataname){
  stopifnot(length(unique(groupname))==1,
            length(unique(dataname))==1,
            length(unique(h5filename))==1)

  fh <- rhdf5::H5Fopen(h5filename,flags = "H5F_ACC_RDWR")
  if(substr(groupname,1,1)!="/"){
    groupname <- paste0("/",groupname)
  }
  groupname <- ifelse(groupname=="/","",groupname)
  d_path <- paste0(groupname,"/",dataname)
  dsp <-rhdf5::H5Oopen(fh,d_path)
  attr <- rhdf5::H5Aopen(dsp,"doTranspose")
  doTranspose <- rhdf5::H5Aread(attr)==0
  rhdf5::H5Aclose(attr)
  rhdf5::H5Oclose(dsp)
  rhdf5::H5Fclose(fh)
  if(doTranspose){
    return(t(rhdf5::h5read(h5filename,d_path)))
  }else{
    return(rhdf5::h5read(h5filename,d_path))
  }
}


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

list.datasets <- function(h5filename,groupname="/",subcols=NULL){
  if(is.null(groupname)){
    groupname <- "/"
  }
  if(substr(groupname,1,1)!="/"){
    groupname <- paste0("/",groupname)
  }
  if(is.null(subcols)){
    return(rhdf5::h5ls(h5filename) %>% dplyr::filter(group==groupname) %>% dplyr::select(name) %>% dplyr::pull(1))
  }else{
    return(rhdf5::h5ls(h5filename) %>% dplyr::filter(group==groupname,name %in% subcols) %>% dplyr::select(name) %>% dplyr::pull(1))
  }
}






read_df_h5 <- function(h5filepath,groupname=NULL,subcols=NULL,filtervec=NULL){
  library(rhdf5)
  stopifnot(file.exists(h5filepath))
  if(is.null(groupname)){
    groupname <- "/"
  }
  if(substr(groupname,1,1)!="/"){
    groupname <- paste0("/",groupname)
  }
  stopifnot(group_exists(h5filepath,groupname))

  dsets <- list.datasets(h5filepath,groupname,subcols)
  if(substr(groupname,nchar(groupname),nchar(groupname))!="/"){
    groupname <- paste0(groupname,"/")
  }
  paths <- as.list(paste0(groupname,dsets))
  names(paths) <- dsets
  stopifnot(length(dsets)>0)
  if(is.logical(filtervec)){
    filtervec <- which(filtervec)
  }
  return(purrr::map_dfc(.x = paths,purrr::compose(c,rhdf5::h5read),file=h5filepath,index=list(filtervec)))
}


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




#' write_mat_h5
#' write a matrix to an HDF5 file
#' @param h5file the name of the HDF5 file to write
#' @param groupname the name of the group for the HDF5 file (can specify several groups, e.g 'a/b/c')
#' @param dataname the name of the matrix in the file
#' @param data the matrix to write, can be of any atomic type
#' @param deflate_level integer specifying compression level
#' @param doTranspose Bool indicating whether or not to write the data in row-major or column major order
write_mat_h5 <- function(h5file, groupname="/", dataname, data, deflate_level = as.integer(c(0)),doTranspose=T){
    library(rhdf5)
    prep_h5file(h5file,create_dir=T)


    groupname <- ifelse(groupname=="/","",groupname)
    if(!group_exists(h5file,groupname)){
        rhdf5::h5createGroup(h5file,groupname)
    }
    d_path <- paste0(groupname,"/",dataname)
    if(doTranspose){
        rhdf5::h5write(data,h5file,d_path)
    }else{
        rhdf5::h5write(t(data),h5file,d_path)
    }

    fh <- rhdf5::H5Fopen(h5file,flags = "H5F_ACC_RDWR")
    dsp <-rhdf5::H5Oopen(fh,d_path)
    transpose_wr <- ifelse(doTranspose,1L,0L)
    rhdf5::h5writeAttribute(transpose_wr,dsp,"doTranspose")
    rhdf5::H5Oclose(dsp)
    rhdf5::H5Fclose(fh)
}

#' write_df_h5
#' Write a dataframe to an HDF5 file
#' @param df dataframe to write (current support for nested/list dataframes is iffy)
#' @param groupname name of group in which to write the dataframe, this will basically be the name of the dataframe in the HDF5 file
#' @param outfile the name of the hdf5 file to write
#' @param deflate_level integer specifying level of compression to apply, with higher indicating higher compression
write_df_h5 <- function(df,groupname="/",outfile,deflate_level=4L){

  prep_h5file(outfile,create_dir = T)


  # dataname <- colnames(df)
  # datapaths <-paste0(groupname,"/",dataname)

  purrr::iwalk(df,
               function(val,name,filename,groupname,deflate_level){
                 write_vec(h5filename = filename,groupname = groupname,dataname = name,data = val,deflate_level=deflate_level)
               },
               filename=outfile,
               groupname=groupname,
               deflate_level=deflate_level)
}


gds2hdf5 <- function(gdsfile,hdf5file,deflate_level=4L){
  library(rhdf5)
  gds <- SeqArray::seqOpen(gdsfile)
  snp_info <-read_SNPinfo_gds(gds,alleles=T,MAF=F,region_id=T,map = T)
  write_df_h5(df = snp_info,
              groupname = "SNPinfo",
              outfile=hdf5file,
              deflate_level = deflate_level)
  dosage2hdf5(gds=gds,hdf5file=hdf5file,deflate_level=deflate_level)



}


dosage2hdf5 <- function(gds,hdf5file,chunksize=5000L,deflate_level=4L){
  library(rhdf5)
  p <- calc_p(gds)
  N <- calc_N(gds)
  dims <- c(p,N)
  is_haplo <- is_haplo(gds)
  if(!dir.exists(dirname(hdf5file))){
    dir.create(dirname(hdf5file))
  }
  if(!file.exists(hdf5file)){
    h5createFile(hdf5file)
  }
  h5createDataset(file = hdf5file,
                  dataset = "dosage",
                  dims = dims,maxdims = dims,
                  storage.mode ="double",chunk = c(chunksize,N),
                  level = deflate_level,
                  showWarnings = F)

  fh <- rhdf5::H5Fopen(hdf5file,flags = "H5F_ACC_RDWR")
  dsp <-rhdf5::H5Oopen(fh,"dosage")
  transpose_wr <- 1L
  rhdf5::h5writeAttribute(transpose_wr,dsp,"doTranspose")
  rhdf5::H5Oclose(dsp)
  rhdf5::H5Fclose(fh)



  write_chunk <-function(index,x,h5loc,is_haplo){
    if(is_haplo){
      tobj <- t(2.0-x)
    }else{
      tobj <- t(x)
    }
    h5write(obj=tobj,
            file=h5loc,
            name="dosage",
            start=c(index,1))
  }





  seqBlockApply(gdsfile = gds,
                var.name = "$dosage",
                FUN =write_chunk,
                h5loc=hdf5file,
                is_haplo=is_haplo,
                var.index="relative")



}
