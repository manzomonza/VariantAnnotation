# AnnotationModules
## Depends on
## -- analysis_dir
## -- snvt

horak_table = "~/github_app/VariantAnnotationModules/HorakScoreTable.txt"
ONCOGEN_POSITIONS = readr::read_tsv("/mnt/NGS_Diagnostik/Variant_databases/OncoKB/Oncokb_Clinvar_oncogenic_positions.tsv")
ONCOGEN_POSITIONS$ProteinChangePosition = stringr::str_extract(string = ONCOGEN_POSITIONS$ProteinChange, pattern = "(?<=\\D)\\d+")


oncogenpositions = table_retrieve_oncogenic_information(parsed_fp$parsed_snv, oncogen_vars = ONCOGEN_POSITIONS)
readr::write_tsv(oncogenpositions, file=annotation_fp$oncogenPos)

snvt = readr::read_tsv(parsed_fp$parsed_snv)
annotation_fp = VariantAnnotationModules::annotation_filepaths(analysis_dir)
write_Annotation_Modules(snvt, annotation_fp = annotation_fp)

### MP Variant Table
source("/home/ionadmin/github/GDrive_VariantReport/Gauths.R")
idoi = googledrive::as_id("https://docs.google.com/spreadsheets/d/1B-NfpRNhadl7w1f5UPkRA_XEg4YI3N4pHRxd9yZgZkc/edit?usp=drive_web&ouid=116704210424700639172")
MPvars = googlesheets4::read_sheet(idoi, skip = 1)
mpvs = MP_check_retrieve_table(snv, MPvars)
mpv_filepath = paste0(analysis_dir,'/annotation_output/annotation_MP_variant.tsv')
if(nrow(mpvs) > 0){
  readr::write_tsv(mpvs, file = mpv_filepath )
}


### BRAF Variant Class
idoi = googledrive::as_id("https://docs.google.com/spreadsheets/d/1xQ3FfHV2JLndT7J_yLHGcrb7Uqpc-cjjio9OYXkjz14/edit")
MPvars = googlesheets4::read_sheet(idoi, skip = 1)
mpvs = BRAF_class_retrieve_table(snv, MPvars)
braf_class = paste0(analysis_dir,'/annotation_output/annotation_BRAF_variant_class.tsv')
if(nrow(mpvs) > 0){
  readr::write_tsv(mpvs, file = braf_class)
}

### HORAK sccores
(horak_scores = Horak_score_function_calls(annotation_fp))
Horak_vals = HorakScore(horak_scores, horakScore_fp = annotation_fp$HorakScoreListings)
Horak_vals$classification = sapply(Horak_vals$Horak_score, Horak_classification)
readr::write_tsv(Horak_vals, file=annotation_fp$Horak)



