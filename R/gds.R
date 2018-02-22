#' Obtain info about variants in a gds object
#'
#'
#' @param gds an `SeqArray` `gds` object (obtained by using the `seqOpen` function in the package `SeqArray`)
#' @param alleles A boolean indicating whether alleles should be returned
#' @param MAF A boolean indicating whether MAF should be computed and returned of not
#' @param region_id A boolean indicating whether to return the LD block of each SNP
read_SNPinfo_gds <- function(gds,alleles=F,MAF=F,region_id=F,map=F,info=F,more=NULL){

  tdf <- tibble::data_frame(SNP=SeqArray::seqGetData(gds,var.name="annotation/id"),
                            snp_id=SeqArray::seqGetData(gds,var.name="variant.id"),
                            chr=SeqArray::seqGetData(gds,var.name="chromosome"),
                            pos=SeqArray::seqGetData(gds,var.name="position"))
  if(alleles){
    tdf <- dplyr::mutate(tdf,allele=SeqArray::seqGetData(gds,var.name="allele"))
  }
  if(MAF){
    library(SeqArray)
    tdf <- dplyr::mutate(tdf,MAF=SeqArray::seqAlleleFreq(gds))
  }
  if(region_id){
    tdf <- dplyr::mutate(tdf,region_id=SeqArray::seqGetData(gds,"annotation/info/LD_chunk"))
  }
  if(map){
    tdf <-dplyr::mutate(tdf,map=SeqArray::seqGetData(gds,"annotation/info/map"))
  }
  if(info){
    tdf <-dplyr::mutate(tdf,info=SeqArray::seqGetData(gds,"annotation/qual")) %>% dplyr::mutate(info=info/max(info))
  }
  if(!is.null(more)){
    tdf <- dplyr::bind_cols(tdf,purrr::map_dfc(more,SeqArray::seqGetData,gdsfile=gds))
  }
  return(tdf)
}


subset_gds <- function(gds,info_df=NULL,region_id=F,sample.id=NULL){
  if(is.null(info_df)){
    si_df <-read_SNPinfo_gds(gds,region_id=region_id) %>% dplyr::distinct(snp_id, .keep_all = T) %>% dplyr::arrange(snp_id)
  }else{
    si_df <- dplyr::inner_join(info_df, read_SNPinfo_gds(gds,region_id=region_id)) %>%
      dplyr::distinct(snp_id, .keep_all = T) %>% dplyr::arrange(snp_id)
  }
  SeqArray::seqSetFilter(gds,variant.id = si_df$snp_id,sample.id=sample.id)
  return(si_df)
}

subset_gds_ind <- function(gds,sample.id){
  seqArray::seqSetFilter(gds,sample.id=sample.id)
}

read_SNPinfo_ldsc_ld <- function(gds){
  return(tibble::data_frame(CHR=SeqArray::seqGetData(gds,var.name="chromosome"),
                            SNP=SeqArray::seqGetData(gds,var.name="annotation/id"),
                            BP=SeqArray::seqGetData(gds,var.name="position"),
                            CM=SeqArray::seqGetData(gds,var.name="annotation/info/map"),
                            MAF=SeqArray::seqAlleleFreq(gds)))
}



read_SNPinfo_ldsc <- function(gds){

  tibble::data_frame(SNP=SeqArray::seqGetData(gds,var.name="annotation/id"),
                            allele=SeqArray::seqGetData(gds,var.name="allele")) %>%
    tidyr::separate(allele,c("A1","A2"))
}


is_SNV <- function(allelevec){
  stopifnot(is.character(allelevec))
  return(sapply(strsplit(allelevec,split=",",fixed=T,useBytes = T),function(x){max(sapply(x,nchar))})==1)
}


subset_export_gds <- function(gds,sample.id=NULL,input_df=NULL,outfile_gds,outfile_h5,reset=T){
  stopifnot(!is.null(outfile_h5),!is.null(outfile_gds),!is.null(sample.id)||is.null(input_df))
  output_df <- subset_gds(gds,input_df,region_id=T,sample.id=sample.id)
  p <- calc_p(gds)
  stopifnot(nrow(output_df)==p)
  gds2hdf5(gds,outfile_h5)
  seqExport(gds,outfile_gds)
  seqClose(gds)
  return(output_df)
}


snpgdsR2SNP <- function(X,snp_df,samp_df=data_frame(sample_id=as.character(1:ncol(X))),outf=tempfile(),compress.geno="LZ4_RA.fast",compress.annotation="LZ4_RA.fast"){
SNPRelate::snpgdsCreateGeno(gds.fn = outf,
                            genmat=X,
                            sample.id = samp_df$sample_id,
                            snp.id = snp_df$snp_id,
                            snp.rs.id =snp_df$SNP,
                            snp.chromosome = snp_df$chr,
                            snp.position = snp_df$pos,
                            snp.allele = snp_df$allele,
                            snpfirstdim = (nrow(X)==nrow(snp_df)),
                            compress.annotation = compress.annotation,
                            compress.geno = compress.geno)
  return(outf)
}


seqIMPUTE2GDS <- function(hap.fn,leg.fn,sample.fn,map.fn,out.gdsfn,chrom,compress.geno="LZ4_RA.fast",compress.annotation="LZ4_RA.fast"){

  map_dat <- readr::read_delim(map.fn,delim=" ",col_names = c("ID","pos","map")) %>% dplyr::mutate(chrom=chrom)
  sample.id <- scan(sample.fn,what=character(),sep="\n")
  sample.id <- c(sapply(sample.id,function(x)c(paste0(x,"-1"),paste0(x,"-2"))))
  leg_df <- readr::read_delim(leg.fn,delim=" ") %>% dplyr::mutate(allele=paste0(allele0,",",allele1),chrom=chrom)
  geno_mat <- readr::read_delim(hap.fn,delim=" ",col_names = sample.id) %>% dplyr::mutate(ID=leg_df$ID)
  leg_df <- dplyr::inner_join(leg_df,map_dat)
  geno_mat <- dplyr::semi_join(geno_mat,leg_df)
  geno_mat <- data.matrix(dplyr::select(geno_mat,-ID))
  p <- nrow(geno_mat)
  tfile <- tempfile()
  SNPRelate::snpgdsCreateGeno(gds.fn = tfile,genmat = geno_mat,
                              sample.id = sample.id,
                              snp.id = 1:p,
                              snp.rs.id = leg_df$ID,
                              snp.chromosome = leg_df$chrom,
                              snp.position = leg_df$pos,
                              snp.allele = leg_df$allele,
                              compress.annotation = compress.annotation,
                              compress.geno = compress.geno)

  SeqArray::seqSNP2GDS(gds.fn =tfile ,out.fn =out.gdsfn ,storage.option ="LZ4_RA.fast")
  gds <- SeqArray::seqOpen(out.gdsfn,readonly = F)
  gdsfmt::add.gdsn(gdsfmt::index.gdsn(gds,"annotation/info"),"map",leg_df$map,replace=T,compress="LZ4_RA.fast")
  SeqArray::seqClose(gds)

}


import_panel_data <- function(temp_gds=NULL,
                              map_df=NULL,
                              output_file=NULL,
                              overwrite=F,
                              parallel=1,
                              ld_break_file=NULL){

  stopifnot(!is.null(output_file),
            !is.null(map_df),
            all(file.exists(temp_gds)),
            !is.null(ld_break_file))

  stopifnot(file.exists(ld_break_file))

  if(file.exists(output_file)){
    stopifnot(overwrite)
    file.remove(output_file)
  }
  SeqArray::seqMerge(temp_gds,output_file,storage.option="LZ4_RA.fast")

  gds <-SeqArray::seqOpen(output_file,readonly = F)
  leg_df <- read_SNPinfo_gds(gds) %>% dplyr::left_join(map_df)
  stopifnot(nrow(leg_df)==length(SeqArray::seqGetData(gds,"variant.id")))


  gdsfmt::add.gdsn(gdsfmt::index.gdsn(gds,"annotation/info"),"map",leg_df$map,replace=T,compress="LZ4_RA.fast")
  cat("Adding Chunk Delimiters\n")
  SeqArray::seqClose(gds)
  add_chunk_gds(output_file,ld_break_file)
}


add_chunk_gds <- function(gds_file,region_bed_file){

  stopifnot(file.exists(region_bed_file),file.exists(gds_file))
  gds <- SeqArray::seqOpen(gds_file,readonly = F)
  region_bed <- readr::read_delim(region_bed_file,delim="\t",trim_ws = T) %>% dplyr::mutate(range_id=as.character(1:n()),chr=gsub("chr","",chr))
  #  region_range <- GenomicRanges::GRanges(seqnames=region_bed$chr,ranges=IRanges::IRanges(start=region_bed$start,end = region_bed$stop))
  snp_df <- read_SNPinfo_gds(gds) %>% dplyr::mutate(snp_id=as.character(snp_id))
  match_df <- match_SNP(region_bed,snp_df) %>% dplyr::arrange(as.integer(snp_id))
  stopifnot(all(snp_df$snp_id==match_df$snp_id))
  gdsfmt::add.gdsn(gdsfmt::index.gdsn(gds,"annotation/info"),"LD_chunk",as.integer(match_df$range_id),replace=T,compress="LZ4_RA.fast")
  SeqArray::seqClose(gds)
}

#' match_SNP
#' Given a dataframe of ranges, match another dataframe of SNPs to that range
#' @param break_df A dataframe having the following columns:
#' chr (character or integer)
#' start (integer)
#' stop (integer)
#' range_id (integer or character)
#' @param snp_df A dataframe having the following columns
#' chr (character or integer)
#' pos (integer)
#' snp_id (character or integer)
#' @param match_at_start indicate whether you want to boundary SNPs to match at start or stop
#' @return
#' A dataframe joining the two dataframes where start<=pos<=stop
match_SNP <- function(break_df,snp_df,match_at_start=T){

  # Whether or not break_df and snp_df use character or integer for chr, we'll use character
  break_cols <- colnames(break_df)
  snp_cols <- colnames(snp_df)
  stopifnot(all(c("chr","start","stop","range_id") %in% break_cols))
  stopifnot(all(c("chr","pos","snp_id") %in% snp_cols))

  snp_chr <- is.character(snp_df$chr)
  break_chr <- is.character(break_df$chr)
  stopifnot((!snp_chr & !break_chr) | (snp_chr & break_chr))

  stopifnot(nrow(dplyr::anti_join(snp_df, break_df, by = "chr")) == 0)
  if (match_at_start) {
    break_df <- dplyr::mutate(break_df,tstop = stop - 1)
    break_df <- dplyr::mutate(break_df,tstart = start)
  }else{
    break_df <- dplyr::mutate(break_df,tstart = start + 1)
    break_df <- dplyr::mutate(break_df,tstop = stop)
  }
  #Start by just subsetting break_df by chr in snp_df
  ch_break_df <- dplyr::semi_join(break_df, snp_df, by = "chr")
  ch_breakl <- split(ch_break_df,ch_break_df$chr)
  snp_breakl <- split(snp_df,snp_df$chr)

  stopifnot(length(ch_breakl) == length(snp_breakl))
  stopifnot(all(names(ch_breakl) == names(snp_breakl)))

    match_df <- dplyr::bind_rows(
                           purrr::map2(
                                      ch_breakl,snp_breakl,
                                      function(ch_break_df,snp_df){

                                          ld_ranges <- IRanges::IRanges(
                                                                    ch_break_df$tstart,
                                                                    ch_break_df$tstop,names = ch_break_df$range_id
                                                                )
                                          pos_ranges <- IRanges::IRanges(snp_df$pos,snp_df$pos,names = snp_df$snp_id)
                                          ol <-   IRanges::findOverlapPairs(pos_ranges,ld_ranges)
                                          match_df <- data_frame(snp_id = ol@first@NAMES,range_id = ol@second@NAMES) %>%
                                              inner_join(snp_df,by = "snp_id") %>%
                                              inner_join(ch_break_df,by = c("range_id","chr"))
                                          return(match_df)
  }))
  n_snp_matches <- group_by(match_df,chr) %>% summarise(nmatches = n())
  n_snp_ct <- group_by(snp_df,chr) %>% summarise(nsnp=n())
  snp_match_sum <- inner_join(n_snp_matches,n_snp_ct,by="chr")
  n_unmapped <- filter(snp_match_sum,nmatches<nsnp) %>% summarise(n_unmapped=sum(nmatches)-sum(nsnp)) %>% pull(n_unmapped)
  n_dups <- filter(snp_match_sum,nmatches>nsnp) %>% summarise(n_dups=sum(nmatches)-sum(nsnp)) %>% pull(n_dups)
  if (n_unmapped > 0){
    warning(paste0(n_unmapped, " SNPs did not map to any region"))
  }
  if (n_dups > 0) {
    warning(paste0(n_dups, " SNPs mapped to multiple regions"))
  }
  match_df <- group_by(match_df,chr) %>% summarise(p=n_distinct(snp_id)) %>% inner_join(match_df)

  return(match_df %>% select(-tstart,-tstop))
}





#' Find number of SNPs in the dataset
#' @param gds A SeqArray gds object
calc_p <- function(gds){
  length(SeqArray::seqGetData(gds,"variant.id"))
}

#' Find number of individuals in the dataset
#' @param gds A SeqArray gds object
calc_N <- function(gds){
  length(SeqArray::seqGetData(gds,"sample.id"))
}


is_haplo <- function(gds,checkSNPs=F){
  samples <- SeqArray::seqGetData(gds,"sample.id")
  # sample_trim <- substr(samples,1,nchar(samples)-2)

  sample_tail <- purrr::map_chr(strsplit(samples,split="-",fixed=T,useBytes=T),purrr::possibly(function(x)x[2],NA_character_,quiet=F))
  if(all(sample_tail %in% c("1","2"))){
    return(T)
  }
  return(F)
}

seqGetHaploDosage <- function(gds){
  stopifnot(is_haplo(gds))
  H <- 2.0-SeqArray::seqGetData(gds,var.name="$dosage")
  stopifnot(max(H)==1)
  return(H)
}



calc_LD_gds <- function(gds,m=85,Ne=11490.672741,cutoff=1e-3){
  map_dat <- SeqArray::seqGetData(gds,var.name="annotation/info/map")
  H <- seqGetHaploDosage(gds)
  return(LDshrink::calcLD(hmata = H,mapa = map_dat,m = m,Ne = Ne,cutoff = cutoff))
}


read_region_id <- function(gdsfile,region_id,center=F){
  gds <- SeqArray::seqOpen(gdsfile)
  filter_region_id(gds,as.integer(region_id))
  dat <-scale(seqGetData(gds,"$dosage"),center=center,scale=F)
  SeqArray::seqClose(gds)
  return(dat)
}


filter_region_id <- function(gds,region_id){
  LD_chunks <- SeqArray::seqGetData(gds,var.name="annotation/info/LD_chunk")
  stopifnot(sum(LD_chunks%in%region_id)>0)
  SeqArray::seqSetFilter(gds,variant.sel=LD_chunks %in% region_id)
}

chunkwise_LDshrink <- function(gds_file,region_id=1,outfile=NULL,m=85,Ne=11490.672741,cutoff=1e-3,evd=T){
  library(rhdf5)
  stopifnot(file.exists(gds_file),!is.null(outfile))
  gds <- SeqArray::seqOpen(gds_file,readonly = T)
  filter_region_id(gds,region_id)

  si <- read_SNPinfo_gds(gds) %>% dplyr::mutate(region_id=region_id)
  R <- calc_LD_gds(gds,m=m,Ne=Ne,cutoff=cutoff)


  write_df_h5(df=si,groupname = "LDinfo",outfile = outfile,deflate_level = 4)
  tls <- h5ls(outfile)
  tls <- paste0(tls$group,"/",tls$name)
  stopifnot(any("/LDinfo/SNP" %in% tls))
  write_mat_h5(outfile,groupname = "LD",dataname = "R",data = R,deflate_level = 4,doTranspose = F)
  if(evd){
    eigenR <- eigen(R)
    stopifnot(min(eigenR$values)>0)
    write_mat_h5(outfile,groupname="EVD",dataname = "Q",data = eigenR$vectors,deflate_level = 0L)
    write_vec(outfile,groupname="EVD",dataname = "D",data = eigenR$values,deflate_level = 2L)
  }
  return(dim(R))
}


chunkwise_LDshrink_ldsc <- function(gds_file,chrom=19,out_dir=NULL,m=85,Ne=11490.672741,cutoff=1e-3){

  stopifnot(!is.null(gds_file),
            file.exists(gds_file),
            !is.null(out_dir),dir.exists(out_dir))
  outfile <- file.path(out_dir,paste0(chrom,".l2.ldscore.gz"))
  soutfile <- file.path(out_dir,paste0(chrom,".l2.M_5_50"))
  gds <-SeqArray::seqOpen(gds_file,readonly = T)
  # chroms <-SeqArray::seqGetData(gds,"chromosome")
  SeqArray::seqSetFilterChrom(gds,include=as.character(chrom))
  LD_chunks <-SeqArray::seqGetData(gds,var.name="annotation/info/LD_chunk")
  region_ids <- unique(LD_chunks)
  mfilt <-SeqArray::seqGetFilter(gds)
  resl <- list()
  # pb <- progress::progress_bar$new(total=length(region_ids))
  for(i in 1:length(region_ids)){
    region_id <- region_ids[i]
    tr <- LD_chunks %in% region_id
    SeqArray::seqSetFilter(gds,variant.sel=tr,action="push+intersect")
    si <- read_SNPinfo_ldsc_ld(gds)
    R <- calc_LD_gds(gds,m = m,Ne = Ne,cutoff = cutoff)
    si <- dplyr::mutate(si,L2=colSums(R^2)-1)
    resl[[i]] <- si
    SeqArray::seqSetFilter(gds,action="pop")
    tfilt <-SeqArray::seqGetFilter(gds)
    stopifnot(all.equal(mfilt,tfilt))
    # pb$tick()
  }
  resdf <- dplyr::bind_rows(resl)
  readr::write_delim(resdf,path=outfile,delim = "\t")
  nc <- dplyr::filter(resdf,MAF>0.05) %>% nrow
  write(x = nc,file = soutfile)
  return(dim(R))
}






