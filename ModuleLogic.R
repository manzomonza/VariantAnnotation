# Module logic

CLINVAR = data.table::fread('/Users/manzo/USB/USB_Diagnostics/VariantAnnotationModules/testing/dbs/snippet_clinvar.txt')
CLINVAR = dplyr::filter(CLINVAR, Assembly == "GRCh37")

COSMIC = data.table::fread('/Users/manzo/USB/USB_Diagnostics/VariantAnnotationModules/testing/dbs/head_cosmic.tsv')
table(COSMIC$`Gene name`)
testfiles = list.files(path = "./testing", recursive = TRUE, full.names = TRUE)
## SNVS

snvfiles = grep("snv.txt", testfiles, value = TRUE)
snvs = lapply(snvfiles, readr::read_tsv)


table(unlist(lapply(snvs[[4]]$protein, oneORthree_code)))

clinvar_check_gene = function(genestr, clinvar_fil){
  clinvar_fil = dplyr::filter(clinvar_fil, GeneSymbol == genestr )
  return(clinvar_fil)
}
clinvar_check_coding = function(codingstr, clinvar_fil){
  clinvar_fil = dplyr::filter(clinvar_fil, grepl(codingstr, Name, fixed = TRUE ))
  return(clinvar_fil)
}

clinvar_check_protein = function(proteinstr, clinvar_fil){
  clinvar_fil = dplyr::filter(clinvar_fil, grepl(proteinstr, Name, fixed = TRUE ))
  return(clinvar_fil)
}



clinvar = clinvar_check_coding(codingstr = codingstr, clinvar_fil = clinvar)

dplyr::filter(CLINVAR, GeneSymbol == "CFTR")

chit = clinvar_filtering(clinvar = CLINVAR,
                     genestr = 'CFTR',
                     proteinstr = 'p.Gln493Ter',
                     codingstr = NA)
clinvar_significance(chit)
clinvar_variantID(chit)

dplyr::count(CLINVAR,GeneSymbol, sort = TRUE)




amino_acid_conversion_three_to_one()
stringr::str_extract(CLINVAR$Name[1:11], pattern = "p\\.\\w+")
stringr::str_extract(CLINVAR$Name[1:11], pattern = "p\\.\\w+")

clinvar_extract_coding(CLINVAR$Name[11:20])
clinvar_extract_transcript(CLINVAR$Name[11:20])
clinvar_extract_protein(CLINVAR$Name[11:20])

colnames(CLINVAR)
CLINVAR$VariationID

grepl('c.892+48G>A', CLINVAR$Name[11:20], fixed = TRUE)
grepl('p.Gly93Arg', CLINVAR$Name[11:20], fixed = TRUE)
grepl('c.25-2A>G', CLINVAR$Name[11:20], fixed = TRUE)


sapply(snvs[[3]]$protein, fsClinvarfix)
lapply(snvs, function(x) x$protein)

snvs = lapply(snvs, amino_acid_code_three_to_one)
## 1. sheet -- overview
lapply(snvs, variants)

## 2. sheet -- technical classification
lapply(snvs, tech_class)


