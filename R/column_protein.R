# ## 0. Adjust protein three letter to one letter
#' Convert snv table with protein column from three letter to single letter amino acid code
#'
#' @param snv
#'
#' @return
#' @export
#'
#' @examples
amino_acid_code_3_to_1 = function(snv){
  snv$protein = unname(sapply(snv$protein, amino_acid_conversion_three_to_one))
  return(snv)
}

#' Convert snv table with protein column from one letter to three letter amino acid code
#'
#' @param snv
#'
#' @return
#' @export
#'
#' @examples
amino_acid_code_1_to_3 = function(snv){
  snv$protein = unname(sapply(snv$protein, amino_acid_conversion_one_to_three))
  return(snv)
}
