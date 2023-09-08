## oncogen_pos 2019

#' Retrieve if variant has been identified as oncogenic position
#'
#' @param gene
#' @param coding_change
#' @param protein_change
#' @param oncogen_pos
#'
#' @return
#' @export
#'
#' @examples
oncogen_pos_retrieve_interpretation = function(gene, coding_change, protein_change, oncogen_pos){
  if(is.na(protein_change)){
    onco_hit = dplyr::filter(oncogen_pos, gene == gene & coding == coding_change)
    if(nrow(onco_hit) >0){
      onco_intp = paste0(unique(onco_hit$oncogenic), collapse = ";")
    }else{
      onco_intp = NA
      }
  }else{
    onco_hit = dplyr::filter(oncogen_pos, gene == gene & ProteinChange == protein_change)
    if(nrow(onco_hit) >0){
      onco_intp = paste0(unique(onco_hit$oncogenic), collapse = ";")
    }else{
      onco_intp = NA
    }
  }
  return(onco_intp)
}

#' Apply oncogenic position check to whole table.
#'
#' @param filepath
#' @param oncogen_pos
#'
#' @return
#' @export
#'
#' @examples
table_retrieve_oncogen_positions = function(filepath, oncogen_pos){
  toi = readr::read_tsv(filepath)
  toi = VariantAnnotationModules::amino_acid_code_3_to_1(toi)
  toi$oncogen_pos = NA
  for(i in 1:nrow(toi)){
    toi$oncogen_pos[[i]] = oncogen_pos_retrieve_interpretation(gene = toi$gene[i],
                                                           coding_change = toi$coding[i],
                                                           protein_change = gsub("p\\.", '', toi$protein[i]),
                                                           oncogen_pos = oncogen_pos)
  }
  toi = dplyr::select(toi, rowid, gene, coding, protein, oncogen_pos)
  return(toi)
}



