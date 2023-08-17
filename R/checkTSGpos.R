

TSG_LENGTHS = readRDS("/Users/manzo/USB/USB_Diagnostics/VariantAnnotationModules/testing/dbs/OncoKB_TSG_maxLength.RDS")

#' Check if Gene is TSG and if variant leads to Ter or frameshift within 90% of AA length
#'
#' @param gene
#' @param aa_pos
#' @param TSG_list
#'
#' @return 'likely pathogenic' if criteria are fullfilled, NA otherwise
#' @export
#'
#' @examples
checkTSG <- function(gene, aa_pos, TSG_list){
  if(gene %in% names(TSG_list) & !is.na(aa_pos)){
    if(aa_pos <= round(0.9 * TSG_list[[gene]])){
      interpretation = 'likely pathogenic'
      return(interpretation)
    }else{
      return(NA)
    }
  }else{
    return(NA)
  }
}


#' tsgParseTable Checks SNV contains 'fs' or 'Ter'
#'
#' @param snvtable
#'
#' @return snvtable with tsgInfo column added
#' @export
#'
#' @examples
tsgParseTable <- function(snvtable, TSG_list){
  if(nrow(snvtable) >0){
    asnv = VariantAnnotationModules::amino_acid_code_3_to_1(snvtable)
    asnv$tsgInfo = NA
    asnv$aa_position = NA
    for (i in 1:nrow(asnv)){
      if(grepl("\\*|fs", asnv$protein[i])){
        asnv$aa_position[i] = extract_number_from_alphanumeric_string(asnv$protein[i])
        asnv$tsgInfo[i] = checkTSG(gene = asnv$gene[i],
                                       aa_pos = asnv$aa_position[i],
                                       TSG_list = TSG_list)

      }
    }
    #snvtable <- subset(snvtable, selec=-aa_position)

    return(asnv)
  }
}
