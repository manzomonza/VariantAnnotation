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
oncogen_var_retrieve_interpretation = function(geneid, coding_change, protein_change, oncogen_vars){
  if(is.na(protein_change)){
    onco_hit = dplyr::filter(oncogen_vars, gene == geneid )
    onco_hit = dplyr::filter(onco_hit, coding == coding_change)
    if(nrow(onco_hit) > 0 ){
      onco_intp = paste0(unique(onco_hit$oncogenic), collapse = ";")
    }else{
      onco_intp = NA
      }
  }else{
    onco_hit = dplyr::filter(oncogen_vars, gene == geneid)
    onco_hit = dplyr::filter(onco_hit, ProteinChange == protein_change)
    if(nrow(onco_hit) > 0 ){
      onco_intp = paste0(unique(onco_hit$oncogenic), collapse = ";")
    }else{
      onco_intp = NA
    }
  }
  return(onco_intp)
}

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
oncogenic_Position = function(geneid, ChangePosition, oncogen_vars){
    onco_hit = dplyr::filter(oncogen_vars, gene == geneid)
    onco_hit = dplyr::filter(onco_hit, ProteinChangePosition == ChangePosition)
    if(nrow(onco_hit) >0){
      oncogenic_position = paste0(unique(onco_hit$ProteinChangePosition), collapse = ";")
    }else{
      oncogenic_position = NA
    }
  return(oncogenic_position)
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
table_retrieve_oncogenic_information = function(filepath, oncogen_vars){
  toi = readr::read_tsv(filepath)
  toi$protein = unname(sapply(toi$protein, VariantStringConversions::amino_acid_conversion_three_to_one))
  toi$ChangePosition = stringr::str_extract(string = toi$protein, pattern = "(?<=\\D)\\d+")
  toi$oncogenic_var = NA
  toi$oncogenic_position = NA
  ####
  for(i in 1:nrow(toi)){
    toi$oncogenic_var[[i]] = oncogen_var_retrieve_interpretation(geneid = toi$gene[i],
                                                           coding_change = toi$coding[i],
                                                           protein_change = gsub("p\\.", '', toi$protein[i]),
                                                           oncogen_vars = oncogen_vars)

    toi$oncogenic_position[[i]] = oncogenic_Position(geneid = toi$gene[i],
                                                     ChangePosition = toi$ChangePosition[i],
                                                     oncogen_vars = oncogen_vars)

  }
  toi = dplyr::select(toi, rowid, gene, coding, protein, oncogenic_var, oncogenic_position)
  return(toi)
}



