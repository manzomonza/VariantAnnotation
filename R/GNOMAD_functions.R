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
gnomad_ref_alt_MAF = function(chr, ref, alt, position, gnomad){
  gnomad_hit = dplyr::filter(gnomad, POS == position)
  gnomad_hit = dplyr::collect(gnomad_hit)
  gnomad_hit = dplyr::filter(gnomad_hit, CHR == chr)
  gnomad_hit = dplyr::filter(gnomad_hit, REF == ref)
  gnomad_hit = dplyr::filter(gnomad_hit, ALT == alt)
  if(nrow(gnomad_hit) != 1){
    return(NA)
  }
  return(gnomad_hit$AF)
}

#' Title
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
gnomad_coding_protein_MAF = function(chr, codingstr, proteinstr, position, gnomad){
  gnomad_hit = dplyr::filter(gnomad, POS == position)
  gnomad_hit = dplyr::filter(gnomad_hit, CHR == chr)
  gnomad_hit = dplyr::collect(gnomad_hit)

  if(is.na(is.na(codingstr)) & is.na(proteinstr)){
    return(NA)
  }else  if(is.na(codingstr)){
    gnomad_hit = dplyr::filter(gnomad_hit, protein == proteinstr)
  }else{
    gnomad_hit = dplyr::filter(gnomad_hit, coding == codingstr)
  }
  if(nrow(gnomad_hit) == 1){
    return(gnomad_hit$AF)
  }else{
    return(NA)
  }
 }



#' Assign identified MAF to tsnv table
#'
#' @param snv_table
#'
#' @return
#' @export
#'
#' @examples
table_retrieve_gnomad_MAF = function(snv_table, gnomad){
  snv_table$gnomad_MAF = NA
  snv_table = VariantAnnotationModules::amino_acid_code_1_to_3(snv_table)
  for(i in 1:nrow(snv_table)){
    maf_v = gnomad_ref_alt_MAF(chr = stringr::str_remove(snv_table$seqnames[i], pattern = "chr"),
                       ref = snv_table$ref[i],
                       alt = snv_table$alt[i],
                       position = snv_table$origPos[i],
                       gnomad = gnomad)
    if(is.na(maf_v)){
      maf_v = gnomad_coding_protein_MAF(chr = stringr::str_remove(snv_table$seqnames[i], pattern = "chr"),
                                codingstr = snv_table$coding[i],
                                proteinstr = snv_table$protein[i],
                                position = snv_table$origPos[i],
                                gnomad = gnomad)
    }
    snv_table$gnomad_MAF[i] = maf_v
  }
  return(snv_table)
}


