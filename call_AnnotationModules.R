# AnnotationModules
## Depends on
## -- analysis_dir
## -- snvt
gdrive_modules = "/home/ionadmin/github_app/VariantAnnotationModules/call_gdrive_modules.R"
horak_table = "~/github_app/VariantAnnotationModules/HorakScoreTable.txt"

ONCOGEN_POSITIONS = readr::read_tsv("/mnt/NGS_Diagnostik/Variant_databases/OncoKB/Oncokb_Clinvar_oncogenic_positions.tsv")
ONCOGEN_POSITIONS$ProteinChangePosition = stringr::str_extract(string = ONCOGEN_POSITIONS$ProteinChange, pattern = "(?<=\\D)\\d+")


oncogenpositions = table_retrieve_oncogenic_information(parsed_fp$parsed_snv, oncogen_vars = ONCOGEN_POSITIONS)
readr::write_tsv(oncogenpositions, file=annotation_fp$oncogenPos)

snvt = readr::read_tsv(parsed_fp$parsed_snv)


annotation_dir = paste0(output_path, "/annotation_output")
if(!dir.exists(annotation_dir)){
  dir.create(annotation_dir)
}
path_annot_gnomad = paste0(annotation_dir, '/annotation_gnomad.tsv')
path_annot_cancerHotspot = paste0(annotation_dir, '/annotation_cancerHotspot.tsv')
path_annot_COSMIC = paste0(annotation_dir, '/annotation_COSMIC.tsv')
path_annot_clinvar = paste0(annotation_dir, '/annotation_ClinVar.tsv')
path_annot_Horak = paste0(annotation_dir, '/annotation_horak.tsv')
path_annot_TSG = paste0(annotation_dir, '/annotation_TSG.tsv')
path_annot_oncogenPos = paste0(annotation_dir, '/annotation_oncogenicPositions.tsv')
path_annot_HorakScoreListings = paste0(annotation_dir, '/annotation_HorakScoreListings.tsv')

## WRITE OUT MODULE OUTPUTS
# Cancer Hotspots
cancerhotspot_s = cancerHotspot_info(snvt, cancerHotspots = CANCERHOTSPOTS)
selected_tb = dplyr::select(cancerhotspot_s, rowid, gene, protein, ref_AA, aa_position, mut_AA, mutation_position_count, mutation_count)
readr::write_tsv(selected_tb, file = path_annot_cancerHotspot)

## Cosmic
cosmic_s = COSMIC_function_call(snv_table = snvt, sql_con_tbl = COSMIC)
selected_tb = dplyr::select(cosmic_s, rowid, gene, coding, protein, contains("COSMIC")  )
readr::write_tsv(selected_tb, file = path_annot_COSMIC)

## Gnomad
gnomad_s = table_retrieve_gnomad_MAF(snvt, gnomad = GNOMAD)
selected_tb = dplyr::select(gnomad_s, rowid, gene, ref,alt, coding, protein, contains("gnomad")  )
readr::write_tsv(selected_tb, file = path_annot_gnomad)

## #CLINVAR
clinvar_s = ClinVar_function_call(snvt, clinvar = CLINVAR)
selected_tb = dplyr::select(clinvar_s, rowid, gene, coding, protein, contains("ClinVar"))
readr::write_tsv(selected_tb, file = path_annot_clinvar)

## TSG info
tsg_s = TSG_check_function_call(snvt, TSG_list = TSG_LENGTHS)
selected_tb = dplyr::select(tsg_s, rowid, gene, coding, protein,TSG,canonical_splicesite, aa_position, protein_length )
readr::write_tsv(selected_tb, file = path_annot_TSG)


annotation_fps = list.files(path = annotation_dir)

### HORAK sccores
(horak_scores = Horak_score_function_calls(annotation_fps))
Horak_vals = HorakScore(horak_scores, horakScore_fp = path_annot_HorakScoreListings)
Horak_vals$classification = sapply(Horak_vals$Horak_score, Horak_classification)
readr::write_tsv(Horak_vals, file=path_annot_Horak)

source(gdrive_modules)




