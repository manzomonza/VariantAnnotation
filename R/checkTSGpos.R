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
TSG_length <- function(gene, TSG_list){
  tsg_l = TSG_list[[gene]]
  return(tsg_l)
}


#' Extracts digits following +/- annotation in codingstring
#' Importantly, only extracts first hit
#'
#' @param codingstring
#'
#' @return
#' @export
#'
#' @examples
extract_splitesite = function(codingstring){
  codingpos = unlist(stringr::str_extract_all(codingstring, pattern = "(?<=\\+|-)\\d{1,}"))
  return(codingpos)
}


#' Determine if splice is at +/- 1 or 2 position relative to coding position
#'
#' @param numberlist
#'
#' @return
#' @export
#'
#' @examples
canonical_splicesite = function(numberlist){
  numbervec = unlist(numberlist)
  plusminus_1_or2 = any(numbervec %in% 1:2)
  return(plusminus_1_or2)
}


#' TSG_check_function_call Checks SNV contains 'fs' or 'Ter'
#'
#' @param snvtable
#'
#' @return snvtable with tsgInfo column added
#' @export
#'
#' @examples
TSG_check_function_call <- function(snvtable, TSG_list){
  if(nrow(snvtable) >0){
    asnv = snvtable
    asnv$protein = unname(sapply(asnv$protein, VariantStringConversions::amino_acid_conversion_three_to_one))
    asnv$TSG = NA
    asnv$canonical_splicesite = NA
    asnv$protein_alteration_site = NA
    asnv$aa_position = NA
    asnv$protein_length = NA

    for (i in 1:nrow(asnv)){
      genestring = as.character(asnv$gene[i])
      if(genestring %in% names(TSG_list)){
        asnv$TSG[i] = TRUE
        if(is.na(asnv$protein[i])){
          if(grepl("splice", asnv$location[i])){
            splicesite = extract_splitesite(asnv$coding[i])
            asnv$canonical_splicesite[i] = canonical_splicesite(splicesite)
          }
        }else{
          if(grepl("\\*|fs", asnv$protein[i])){
            asnv$aa_position[i] = extract_number_from_alphanumeric_string(asnv$protein[i])
            asnv$protein_length[i] = TSG_list[[genestring]]
          }else{
            next
          }
        }
      }else{
        next
      }
    }
    return(asnv)
  }
}

