# ## HORAK scoring

#' Apply Horak scoring rule to gnomad MAF
#' if MAF > 5%, -8
#' if 1% < MAF < 5%, -4
#' absent or < 1%, 1
#' @param gnomadpath
#'
#' @return
#' @export
#'
#' @examples
Horak_score_gnomad = function(gnomadpath){
  gnomad = readr::read_tsv(gnomadpath)
  mafs = gnomad$gnomad_MAF
  hscore = ifelse(is.na(mafs), 1,
                  ifelse(mafs > 0.05, -8,
                         ifelse(mafs > 0.01 & mafs < 0.05, -4, 1)))
  gnomad$gnomad_hscore = hscore
  return(gnomad)
}

#lapply(df, function(x) Horak_score_gnomad(x$paths[4]))

#' Apply Horak scoring rule to considering cancerHotspot mutation counts
#'
#' @param oncogenPos_path
#'
#' @return
#' @export
#'
#' @examples
Horak_score_oncogenpos = function(oncogenPos_path){
  oncogenPos_df = readr::read_tsv(oncogenPos_path)
  oncogenPos_df$oncPos_hscore = ifelse(!is.na(oncogenPos_df$oncogenic_var), 4,
                                       ifelse(!is.na(oncogenPos_df$oncogenic_position), 2, 0))

  return(oncogenPos_df)
}

#' Apply Horak scoring rule to considering cancerHotspot mutation counts
#'
#' @param cancerhotspotpath
#'
#' @return
#' @export
#'
#' @examples
Horak_score_TSG = function(TSGpath){
  tsg = readr::read_tsv(TSGpath)
  tsg$TSG_hscore = NA
  for (i in 1:nrow(tsg)){
    if(!is.na(tsg$TSG[i])){
      aa_position = as.numeric(tsg$aa_position[i])
      protein_length = as.numeric(tsg$protein_length[i])
      length_ratio = aa_position/protein_length <= 0.95
      if(any(c(length_ratio, tsg$canonical_splicesite[i]) == TRUE)){
        hscore = 8
      }else{
        hscore = 0
      }
    }else{
      hscore = 0
    }
    tsg$TSG_hscore[i] = hscore
  }
  return(tsg)
}


#' Apply Horak scoring rule to considering cancerHotspot mutation counts
#'
#' @param cancerhotspotpath
#'
#' @return
#' @export
#'
#' @examples
Horak_score_cancerHotspot_counts = function(cancerhotspotpath){
  ch = readr::read_tsv(cancerhotspotpath)
  hscore = ifelse(is.na(ch$mutation_position_count) & is.na(ch$mutation_count),0,
                  ifelse(ch$mutation_position_count < 50 & ch$mutation_count >= 10, 2,
                         ifelse(ch$mutation_position_count >= 50 & ch$mutation_count >= 10, 4,
                                ifelse(ch$mutation_count <= 10, 1, 0))))
  ch$cancerHotspotCount_Hscore = hscore
  return(ch)
}


#' Extract paths for annotation tables and call Horak scoring functions
#'
#' @param annotation_paths
#'
#' @return
#' @export
#'
#' @examples
Horak_score_function_calls = function(annotation_paths){
  gnomad_path = grep("annotation_gnomad.tsv", annotation_paths, value = TRUE)
  cancerhotspot_path = grep("annotation_cancerHotspot.tsv", annotation_paths, value = TRUE)
  tsg_path = grep("annotation_TSG.tsv", annotation_paths, value = TRUE)
  oncogenpos_path = grep("annotation_oncogenicPositions.tsv", annotation_paths, value = TRUE)

  if(!identical(gnomad_path, character(0))){
    gnomad_df = Horak_score_gnomad(gnomad_path)
  }
  if(!identical(cancerhotspot_path, character(0))){
    chc_df = Horak_score_cancerHotspot_counts(cancerhotspot_path)
  }
  if(!identical(tsg_path, character(0))){
    tsg_df = Horak_score_TSG(tsg_path)
  }
  if(!identical(oncogenpos_path, character(0))){
    oncogenpos_df = Horak_score_oncogenpos(oncogenpos_path)
  }

  horak_scores = list(gnomad_df, chc_df, tsg_df, oncogenpos_df)

  return(horak_scores)
}

#' Sum up individual module scores to generate Horak Score
#'
#' @param horak_scores
#'
#' @return
#' @export
#'
#' @examples
HorakScore = function(horak_scores, horakScore_fp){
  horak_scores = lapply(horak_scores, function(x) dplyr::select(x, rowid, gene, contains("Hscore")))
  hscores = purrr::reduce(horak_scores, dplyr::left_join, by = c("rowid", "gene"))
  hscores = dplyr::group_by(hscores, rowid, gene)
  readr::write_tsv(hscores, horakScore_fp)
  hscores = dplyr::mutate(hscores, Horak_score = sum(gnomad_hscore, cancerHotspotCount_Hscore, TSG_hscore, oncPos_hscore, na.rm = TRUE ))
  Horak_score_df = dplyr::select(hscores, rowid, gene, Horak_score )
  #hscores = lapply(hscores, function(x) tidyr::pivot_longer(x, -c(gene,rowid), names_to = "category", values_to = "HorakScore"))
  # hscores = dplyr::bind_rows(hscores)
  # hscores = dplyr::reframe(hscores, rowid, gene, HorakScore = sum(HorakScore, na.rm = TRUE))
  # hscores = dplyr::distinct(hscores)
  return(Horak_score_df)
}

#' Assign Horak Score classifications to values
#'
#' @param HorakScore
#'
#' @return
#' @export
#'
#' @examples
Horak_classification = function(HorakScore){
  if(HorakScore < -7){
    return("benign")
  }else if(HorakScore >= -6 & HorakScore <= -1){
    return("likely benign")
  }else if(HorakScore >= 0 & HorakScore <= 5){
    return("VUS")
  }else if(HorakScore >= 6 & HorakScore <= 9){
    return("likely oncogenic")
  }else{
    return('oncogenic')
  }
}
