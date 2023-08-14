## gnomad functions





#' Extract Minor Allel Frequency from gnomad table
#'
#' @param chr
#' @param ref
#' @param alt
#' @param position
#' @param gnomad
#'
#' @return
#' @export
#'
#' @examples
gnomad_MAF = function(chr, ref, alt, position, gnomad){
  gnomad_hit = dplyr::filter(gnomad, POS == position)
  gnomad_hit = dplyr::filter(gnomad_hit, CHROM == chr)
  gnomad_hit = dplyr::filter(gnomad_hit, REF == ref)
  gnomad_hit = dplyr::filter(gnomad_hit, ALT == alt)
  if(nrow(gnomad_hit) != 1){
    return(NA)
  }
  return(gnomad_hit$AF)
}


#' Assign identified MAF to tsnv table
#'
#' @param snv_table
#'
#' @return
#' @export
#'
#' @examples
attach_gnomad_MAF = function(snv_table){
  snv_table$gnomad_MAF = NA
  for(i in 1:nrow(snv_table)){
  maf_v = gnomad_MAF(chr = stringr::str_remove(snv_table$seqnames[i], pattern = "chr"),
                     ref = snv_table$ref[i],
                     alt = snv_table$alt[i],
                     position = snv_table$origPos[i],
                     gnomad = GNOMAD)
  snv_table$gnomad_MAF[i] = maf_v
  }
  return(snv_table)
}


