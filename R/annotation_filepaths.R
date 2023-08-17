#' Generate list of filepaths from annotation folder
#'
#' @param output_path
#'
#' @return
#' @export
#'
#' @examples
annotation_filepaths = function(output_path){
  annotation_dir = paste0(output_path, "/annotation_output")
  if(!dir.exists(annotation_dir)){
    dir.create(annotation_dir)
  }
  filepaths = list(gnomad = paste0(annotation_dir,'/annotation_gnomad.tsv'),
                   cancerHotspot = paste0(annotation_dir,'/annotation_cancerHotspot.tsv'),
                   COSMIC = paste0(annotation_dir,'/annotation_COSMIC.tsv'),
                   clinvar = paste0(annotation_dir,'/annotation_ClinVar.tsv'),
                   Horak = paste0(annotation_dir,'/annotation_horak.tsv'),
                   TSG = paste0(annotation_dir,'/annotation_TSG.tsv'),
                   filepaths = paste0(annotation_dir,'/annotation_filepaths.tsv'))
  return(filepaths)
}
