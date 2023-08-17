#' Count variant occurences per tissue in COSMIC db
#'
#' @param snv_table
#' @param sql_con_tbl
#'
#' @return
#' @export
#'
#' @examples
COSMIC_function_call = function(snv_table, sql_con_tbl){
  asnv = VariantAnnotationModules::amino_acid_code_3_to_1(snv_table)
  asnv$COSMIC_n_total = NA
  asnv$COSMIC_n_tissue = NA
  for(i in 1:nrow(asnv)){
    gene = asnv$gene[i]
    variants_per_tissue = dplyr::filter(sql_con_tbl, gene_name == gene)
    if(is.na(asnv$protein[i]) | asnv$protein[i] == "p.?"){
      variants_per_tissue = dplyr::filter(variants_per_tissue, mutation_cds == asnv$coding[i])
    }else{
      variants_per_tissue = dplyr::filter(variants_per_tissue, mutation_aa == asnv$protein[i])
    }

    variants_per_tissue = dplyr::collect(variants_per_tissue)
    variants_per_tissue = dplyr::distinct(variants_per_tissue)
    if(nrow(variants_per_tissue) > 0){
      asnv$COSMIC_n_total[i] = nrow(variants_per_tissue)
      count_variants_per_tissue = dplyr::count(variants_per_tissue, primary_site, sort = TRUE)
      asnv$COSMIC_n_tissue[i] = paste0(paste0(count_variants_per_tissue$primary_site, " (", count_variants_per_tissue$n, ")"), collapse = "; ")
    }
  }
  return(asnv)
}



#' Write out cancer Hotspot info table
#'
#' @param CH_snv
#' @param CH_info_file
#'
#' @return
#' @export
#'
#' @examples
write_out_cosmic_info = function(COSMIC_snv, COSMIC_info_file){
  selected_snv = dplyr::select(COSMIC_snv, rowid, gene, coding, protein, contains("COSMIC")  )
  readr::write_tsv(selected_snv, file = COSMIC_info_file)
  return(selected_snv)
}

