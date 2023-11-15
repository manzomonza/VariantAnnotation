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
  asnv = snv_table
  asnv$protein = unname(sapply(asnv$protein, VariantStringConversions::amino_acid_conversion_three_to_one))
  asnv$COSMIC_n_total = NA
  asnv$COSMIC_n_tissue = NA

  unique_gene_names = unique(asnv$gene)
  unique_coding = unique(asnv$coding)
  unique_protein = unique(asnv$protein)

  ## Collect SQL database search to improve speed
  sql_con_tbl =  dplyr::filter(sql_con_tbl, gene_name %in% unique_gene_names)
  sql_con_tbl =  dplyr::filter(sql_con_tbl, mutation_cds %in% unique_coding | mutation_aa %in% unique_protein)
  sql_con_tbl = dplyr::collect(sql_con_tbl)

  for(i in 1:nrow(asnv)){
    gene = asnv$gene[i]
    protein = asnv$protein[i]
    coding = asnv$coding[i]
    variants_per_tissue = dplyr::filter(sql_con_tbl, gene_name == gene)
    if(is.na(asnv$protein[i]) | asnv$protein[i] == "p.?"){
      variants_per_tissue = dplyr::filter(variants_per_tissue, mutation_cds == coding)
    }else{
      variants_per_tissue = dplyr::filter(variants_per_tissue, mutation_aa == protein)
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

