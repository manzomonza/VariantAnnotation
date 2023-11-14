## Annotate CNVs
# Use Onco KB CNV annotation to interpret CNV effects
# Threshold with 4x as gain and 0 - <1 as loss
# then use OncoKB LUT


#' Annotate CNV effect
#'
#' @param cnv_table
#' @param cnv_lut
#'
#' @return
#' @export
#'
#' @examples
annotate_cnv = function(cnv_table, cnv_lut, tsg_list){
  cnv_rows = 1:nrow(cnv_table)
  cnv_table$cnv_status = unlist(lapply(cnv_rows, function(i) determine_cnv_status(perc_05 = cnv_table$perc_0.05[i],
                                                                                  perc_95 = cnv_table$perc_0.95[i])))

  cnv_table$cnv_effect = unlist(lapply(cnv_rows, function(i) gene_cnv_effect(gene = cnv_table$gene[i],
                                                                             cnv_status = cnv_table$cnv_status[i],
                                                                             cnv_lut = cnv_lut)))

  cnv_table$TSG_cnv = unlist(lapply(cnv_rows, function(i) tsg_cnv_effect(gene = cnv_table$gene[i],
                                                                             cnv_status = cnv_table$cnv_status[i],
                                                                         tsg_list = tsg_list)))

  cnv_table$HRR_gene = unlist(lapply(cnv_rows, function(i) hrr_cnv(gene = cnv_table$gene[i],
                                                                          hrrgenes = HRRgenes$gene)))
  return(cnv_table)
}

#' Determine CNV status of CNVs
#'
#' @param perc_05
#' @param perc_95
#'
#' @return
#' @export
#'
#' @examples
determine_cnv_status = function(perc_05, perc_95){
  if(perc_95 < 1){
    return("loss")
  }
  if(perc_05 >= 4){
    return("gain")
  }
  return(NA)
}

#' Compare CNV status with CNV effect in cnv lut
#'
#' @param gene
#' @param cnv_status
#' @param cnv_lut
#'
#' @return
#' @export
#'
#' @examples
gene_cnv_effect = function(gene, cnv_status, cnv_lut){
  cnv_lut = dplyr::filter(cnv_lut, GeneSymbol == gene)
  if(!is.na(cnv_status)){
    if(cnv_status == "gain"){
      cnv_lut = dplyr::filter(cnv_lut, ProteinChange == "Amplification")
    }else if(cnv_status == "loss"){
      cnv_lut = dplyr::filter(cnv_lut, ProteinChange == "Deletion")
    }
    if(nrow(cnv_lut) >0){
      cnv_effect = paste(cnv_lut$oncogenic, collapse = ";")
      return(cnv_effect)
    }else{
      return(NA)
    }
  }else{
    return(NA)
  }
}

#' Determine if CNV loss is happening in TSG
#'
#' @param gene
#' @param cnv_status
#' @param tsg_list
#'
#' @return
#' @export
#'
#' @examples
tsg_cnv_effect = function(gene, cnv_status, tsg_list){
  if(is.na(cnv_status)){
    return(NA)
  }
  tsggene = gene %in% names(tsg_list)
  if(tsggene){
    if(cnv_status == "loss"){
      return("likely pathogenic")
    }else{
      return(NA)
    }
  }else{
    return(NA)
  }
}

#' Indicate if gene is present in HRR gene list
#'
#' @param gene
#' @param hrrgenes
#'
#' @return
#' @export
#'
#' @examples
hrr_cnv = function(gene, hrrgenes){
  hrr_gene = gene %in% hrrgenes
  if(hrr_gene){
    return(TRUE)
  }else{
    return(NA)
  }
}
