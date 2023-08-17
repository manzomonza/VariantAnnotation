# AnnotationModules
## Depends on
## -- analysis_dir
## -- snvt

snvt = readr::read_tsv(parsed_fp$parsed_snv)
annotation_fp = annotation_filepaths(analysis_dir)
Annotation_Modules(snvt, annotation_fp = annotation_fp)

