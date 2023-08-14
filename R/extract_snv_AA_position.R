#' Extract position of SNVs
#'
#' @param clinvar_ready_AA
#'
#' @return SNV amino-acid position
#' @export
#'
#' @examples
extract_number_from_alphanumeric_string <- function(alphanumstring){
  if(is.character(alphanumstring) & !is.na(alphanumstring) & !grepl("_", alphanumstring)){
    numstring = stringr::str_extract(string = alphanumstring, pattern = "(?<=\\D)\\d+")
    numstring = as.integer(numstring)
    return(numstring)
  }else{
    return(NA)

  }
}




