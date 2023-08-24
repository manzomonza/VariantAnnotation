# AnnotationModules
## Depends on
## -- analysis_dir
## -- snvt

snvt = readr::read_tsv(parsed_fp$parsed_snv)
annotation_fp = VariantAnnotationModules::annotation_filepaths(analysis_dir)
write_Annotation_Modules(snvt, annotation_fp = annotation_fp)

# horak_table = "/Users/manzo/USB/USB_Diagnostics/VariantAnnotationModules/HorakScoreTable.txt"
#
# testpath = '/Users/manzo/USB/USB_Diagnostics/ShinyVariants/testfiles'
# df = data.frame(paths = list.files(path = testpath, pattern = "annotation_", recursive = TRUE, full.names = TRUE))
#
# df$sample = basename(dirname(dirname(df$paths)))
# df=tibble::as_tibble(df[,c(2,1)])
# df = dplyr::group_by(df, sample)
# df = dplyr::group_split(df)
# (horak_scores = Horak_score_function_calls(df[[7]]$paths))
# HorakScore(horak_scores)
# Horak_classification()
#
# Horak_vals = HorakScore(horak_scores)
# Horak_vals$classification = sapply(Horak_vals$Horak_score, Horak_classification)
