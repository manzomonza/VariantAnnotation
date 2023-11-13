#### HRR CHECK


#' Check if gene is in HRR gene list
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
HRR_check = function(genestr, mp_vars_table){
  return(genestr %in% mp_vars_table$gene)
}


#' HRR gene check table retrieve
#'
#' @param variant_table
#' @param mp_vars_table
#'
#' @return
#' @export
#'
#' @examples
HRR_check_retrieve_table = function(variant_table, mp_vars_table){
  row_indeces = list()
  for(i in 1:nrow(variant_table)){
    hrr_check = HRR_check(genestr = variant_table$gene[i], mp_vars_table = mp_vars_table)
    if(hrr_check){
      row_indeces[[i]] = variant_table$rowid[i]
    }
  }
  mpv = dplyr::filter(variant_table, rowid %in% unlist(row_indeces) & gene %in% mp_vars_table$gene)
  mpv = dplyr::select(mpv, rowid, gene, coding, protein, AF)
  mpv$category = "HRR gene"
  return(mpv)
}



