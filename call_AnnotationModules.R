# AnnotationModules
## Depends on
## -- analysis_dir
## -- snvt

snvt = readr::read_tsv(parsed_fp$parsed_snv)
annotation_fp = VariantAnnotationModules::annotation_filepaths(analysis_dir)
write_Annotation_Modules(snvt, annotation_fp = annotation_fp)

horak_table = "~/github_app/VariantAnnotationModules/HorakScoreTable.txt"
#
# testpath = '/Users/manzo/USB/USB_Diagnostics/ShinyVariants/testfiles'
# df = data.frame(paths = list.files(path = testpath, pattern = "annotation_", recursive = TRUE, full.names = TRUE))
#
# df$sample = basename(dirname(dirname(df$paths)))
# df=tibble::as_tibble(df[,c(2,1)])
# df = dplyr::group_by(df, sample)
# df = dplyr::group_split(df)
(horak_scores = Horak_score_function_calls(annotation_fp))

Horak_vals = HorakScore(horak_scores)
Horak_vals$classification = sapply(Horak_vals$Horak_score, Horak_classification)

readr::write_tsv(Horak_vals, file=annotation_fp$Horak)
