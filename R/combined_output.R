## Generate combined output file using other files


## 0. Adjust output
adjust_input =function(snvdf){
  # if not genexus, then correct amino acid code to one letter
}

#' 1. Sheet -- Overview of detected variantas
#'
#' @param snvdf
#'
#' @return
#' @export
#'
#' @examples
variants = function(snvdf){
  variantsdf = dplyr::select(snvdf, rowid, transcript, seqnames, gene, coding, protein, AF,  exon, location, variant_type)
  return(variantsdf)
}


#' 2. Sheet -- Technical classification of variants, e.g. if artefacts or not
#'
#' @param snvdf
#'
#' @return
#' @export
#'
#' @examples
tech_class = function(snvdf){
  variantsdf = dplyr::select(snvdf, rowid, seqnames, gene, coding, protein, AF, QUAL, totalDepth,FSAF,  FSAR,  FSRF,  FSRR)
  return(variantsdf)
}

#' 3. Sheet -- Interpretation tools
#'
#' @param snvdf
#'
#' @return
#' @export
#'
#' @examples
tech_class = function(snvdf){
  variantsdf = dplyr::select(snvdf, rowid, seqnames, gene, coding, protein, AF, QUAL, totalDepth,FSAF,  FSAR,  FSRF,  FSRR)
  return(variantsdf)
}




