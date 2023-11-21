#' Retrieve if variant codon position is covered in cancerHotspots v2
#'
#' @param gene
#' @param snv_position
#'
#' @return
#' @export
#'
#' @examples
#' Retrieves No. of total entries in CANCER_HOTSPOTS
#' KRAS p.G12C retrieves n_total == 2175
cancerHotspots_mutation_positionCount <- function(geneName, mutation_position, cancerHotspots){
  if(geneName %in% cancerHotspots$Hugo_Symbol){
    ch = dplyr::filter(cancerHotspots, Hugo_Symbol == geneName)
    if(mutation_position %in% ch$Amino_Acid_Position){
      ch = dplyr::filter(ch, Amino_Acid_Position == mutation_position)
      ntotal = unique(dplyr::pull(ch, Mutation_Count))
    return(as.numeric(ntotal))
    }else{
      return(NA)
    }
  }else{
    return(NA)
  }
}

#' Title
#'
#' @param geneName
#' @param mutation_position
#' @param mut_AA
#' @param cancerHotspots
#'
#' @return
#' @export
#'
#' @examples
cancerHotspots_mutationCount <- function(geneName, mutation_position, mut_AA, cancerHotspots){
  if(geneName %in% cancerHotspots$Hugo_Symbol){
    ch = dplyr::filter(cancerHotspots, Hugo_Symbol == geneName)
    if(mutation_position %in% ch$Amino_Acid_Position){
      ch = dplyr::filter(ch, Amino_Acid_Position == mutation_position)
      ch = dplyr::filter(ch, grepl(mut_AA, Variant_Amino_Acid))
      if(nrow(ch) == 1){
        ntotal = dplyr::pull(ch, Variant_Amino_Acid)
        return(ntotal)
      }else{
        return(NA)
      }
    }else{
      return(NA)
    }
  }else{
    return(NA)
  }
}



#' Assign identified MAF to tsnv table
#'
#' @param snv_table
#'
#' @return
#' @export
#'
#' @examples
cancerHotspot_add_columns_to_snv = function(snv_table){
  snv_table$ref_AA = NA
  snv_table$aa_position = NA
  snv_table$mut_AA = NA
  snv_table$mutation_position_count = NA
  snv_table$mutation_count = NA
  for(i in 1:nrow(snv_table)){
    if(is.character(snv_table$protein[i]) &
       !is.na(snv_table$protein[i]) &
       !grepl("_", snv_table$protein[i]) &
       snv_table$protein[i] != 'p.'){
    snv_table$ref_AA[i] = ref_amino_acid(snv_table$protein[i])
    snv_table$aa_position[i] = extract_number_from_alphanumeric_string(snv_table$protein[i])
    snv_table$mut_AA[i] = mutation_amino_acid(snv_table$protein[i])
    }
  }
  return(snv_table)
}


#' Extract reference amino acid from 'appropriate' variant entries
#'
#' @param aa_change
#'
#' @return
#' @export
#'
#' @examples
ref_amino_acid = function(aa_change){
  if(is.character(aa_change) &
     !is.na(aa_change) &
     !grepl("_", aa_change)){
    aa_change = stringr::str_extract(aa_change, pattern = "(?<=p\\.)\\w")
    return(aa_change)
  }
}

#' Extract mutation amino acid from 'appropriate' variant entries
#'
#' @param aa_change
#'
#' @return
#' @export
#'
#' @examples
mutation_amino_acid = function(aa_change){
  if(is.character(aa_change) & !is.na(aa_change) & !grepl("_", aa_change)){
    aa_change = stringr::str_extract(aa_change, pattern = "(?<=p\\.\\w\\d{1,5})\\D")
    return(aa_change)
  }
}

#' Assign identified MAF to tsnv table
#'
#' @param snv_table
#'
#' @return
#' @export
#'
#' @examples
cancerHotspot_add_mutation_values = function(snv_table, cancerHotspots){
  snv_table$mutation_position_count = NA
  snv_table$mutation_count = NA
  for(i in 1:nrow(snv_table)){
    if(is.character(snv_table$protein[i]) &
       !is.na(snv_table$protein[i]) &
       !grepl("_", snv_table$protein[i]) &
       snv_table$protein[i] != 'p.'){

      snv_table$mutation_position_count[i] = cancerHotspots_mutation_positionCount(geneName = snv_table$gene[i],
                                                                                   mutation_position = snv_table$aa_position[i],
                                                                                   cancerHotspots = cancerHotspots)

      snv_table$mutation_count[i] = cancerHotspots_mutationCount(geneName = snv_table$gene[i],
                                                                 mutation_position = snv_table$aa_position[i],
                                                                 mut_AA = snv_table$mut_AA[i],
                                                                 cancerHotspots = cancerHotspots)
      snv_table$mutation_count[i] = extract_number_from_alphanumeric_string(snv_table$mutation_count[i])
      }
  }
  return(snv_table)
}


#' Run cancerHotspot data extraction pipeline
#'
#' @param snv_table
#' @param cancerHotspots
#'
#' @return
#' @export
#'
#' @examples
cancerHotspot_info = function(snv_table, cancerHotspots){
  asnv = snv_table
  asnv$protein = unname(sapply(asnv$protein, VariantStringConversions::amino_acid_conversion_three_to_one))
  asnv = cancerHotspot_add_columns_to_snv(asnv)
  asnv = cancerHotspot_add_mutation_values(asnv, cancerHotspots = cancerHotspots)
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
write_out_cancerHotspot_info = function(CH_snv, CH_info_file){
  selected_snv = dplyr::select(CH_snv, rowid, gene, protein, ref_AA, aa_position, mut_AA, mutation_position_count, mutation_count)
  readr::write_tsv(selected_snv, file = CH_info_file)
  return(selected_snv)
}

