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
MP_variant_check = function(rowid, genestr, codingstr, protein, mp_vars_table){
  if(!is.na(protein) & !identical(protein, 'p.?')){
    protein = VariantStringConversions::amino_code_conversion_three_to_one(protein)
    mpv = dplyr::filter(mp_vars_table, gene == genestr & amino_acid_1l == protein)

  }else if(!is.na(codingstr)){
    mpv = dplyr::filter(mp_vars_table, gene == genestr & coding == codingstr)
  }
  mpv$rowid = rowid
  mpv = dplyr::relocate(mpv, rowid)
  return(mpv)
}

#' Apply MP_variant_check at table level
#'
#' @param variant_table
#' @param mp_vars_table
#'
#' @return
#' @export
#'
#' @examples
MP_check_retrieve_table = function(variant_table, mp_vars_table){
  variant_hits = list()
  for(i in 1:nrow(variant_table)){
    variant_hits[[i]] = MP_variant_check(rowid = variant_table$rowid[i],
                                         genestr = variant_table$gene[i],
                                         protein = variant_table$protein[i],
                                         codingstr = variant_table$coding[i],
                                         mp_vars_table)
  }
  mpv = dplyr::bind_rows(variant_hits)
  return(mpv)
}
