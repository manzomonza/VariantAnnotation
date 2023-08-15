#


snvt = snvs[[1]]
asnv = adjust_protein(snv_table)

## Cancer Hotspots
cancerhotspot_s = cancerHotspot_info(snvt, cancerHotspots = CANCERHOTSPOTS)
selected_snv = dplyr::select(cancerhotspot_s, rowid, gene, protein, ref_AA, aa_position, mut_AA, mutation_position_count, mutation_count)

## Cosmic
cosmic_s = COSMIC_function_call(snv_table = snvt,sql_con_tbl = COSMIC)
selected_snv = dplyr::select(cosmic_s, rowid, gene, coding, protein, contains("COSMIC")  )

## Gnomad
gnomad_s = attach_gnomad_MAF(snvt, gnomad = GNOMAD)
selected_snv = dplyr::select(gnomad_s, rowid, gene, ref,alt, coding, protein, contains("gnomad")  )

