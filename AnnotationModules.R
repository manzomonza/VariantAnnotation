#
COSMIC
CANCERHOTSPOTS
GNOMAD
CLINVAR

snvt = snvs[[4]]
asnv = amino_acid_code_3_to_1(snvt)
asnv = amino_acid_code_1_to_3(snvt)
(dfme = data.frame(one = asnv$protein, three =  snvt$protein))


clinvar_filtering(genestr = 'AP5Z1', codingstr ='c123', proteinstr = 'fff',
                  clinvar = CLINVAR)


annotation_fp = annotation_filepaths('/Users/manzo/USB/USB_Diagnostics/VariantAnnotationModules/testing/')

## Cancer Hotspots
cancerhotspot_s = cancerHotspot_info(snvt, cancerHotspots = CANCERHOTSPOTS)
selected_snv = dplyr::select(cancerhotspot_s, rowid, gene, protein, ref_AA, aa_position, mut_AA, mutation_position_count, mutation_count)
readr::write_tsv(selected_snv, file = annotation_fp$cancerHotspot)

## Cosmic
cosmic_s = COSMIC_function_call(snv_table = snvt,sql_con_tbl = COSMIC)
selected_snv = dplyr::select(cosmic_s, rowid, gene, coding, protein, contains("COSMIC")  )
readr::write_tsv(selected_snv, file = annotation_fp$COSMIC)

## Gnomad
gnomad_s = attach_gnomad_MAF(snvt, gnomad = GNOMAD)
selected_snv = dplyr::select(gnomad_s, rowid, gene, ref,alt, coding, protein, contains("gnomad")  )
readr::write_tsv(selected_snv, file = annotation_fp$gnomad)

## #CLINVAR
clinvar_s = ClinVar_function_call(snvt, clinvar = CLINVAR)
selected_snv = dplyr::select(clinvar_s, rowid, gene, coding, protein, contains("ClinVar"))
readr::write_tsv(selected_snv, file = annotation_fp$clinvar)

