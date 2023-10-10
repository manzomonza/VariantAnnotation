## Read in and parse MP_variant_table

#' Check MP_variant_table for variant interpretations
#'
#' @param rowid
#' @param genestr
#' @param codingstr
#' @param protein
#' @param mp_vars_table
#'
#' @return
#' @export
#'
#' @examples
BRAF_variant_check = function(rowid, genestr, protein, mp_vars_table){
  if(genestr == "BRAF")
    if(!is.na(protein)){
      protein = VCFparse::amino_acid_conversion_three_to_one(protein)
      mpv = dplyr::filter(mp_vars_table, variant == protein)
      if(nrow(mpv) >0){
        mpv$rowid = rowid
        mpv = dplyr::relocate(mpv, rowid)
        return(mpv)
      }
    }
}


#' BRAF class table retrieve
#'
#' @param variant_table
#' @param mp_vars_table
#'
#' @return
#' @export
#'
#' @examples
BRAF_class_retrieve_table = function(variant_table, mp_vars_table){
  variant_hits = list()
  for(i in 1:nrow(variant_table)){
    variant_hits[[i]] = BRAF_variant_check(rowid = variant_table$rowid[i],
                                         genestr = variant_table$gene[i],
                                         protein = variant_table$protein[i],
                                         mp_vars_table)
  }
  mpv = dplyr::bind_rows(variant_hits)
  return(mpv)
}



