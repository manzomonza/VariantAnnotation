
GNOMAD = data.table::fread('/Users/manzo/USB/USB_Diagnostics/VariantAnnotationModules/testing/dbs/gnomad_10k.tsv')
tsnv = readr::read_tsv('/Users/manzo/USB/USB_Diagnostics/VariantAnnotationModules/testing/test_snv.tsv')

tsnv$seqnames[1:3] = 'chr19'
tsnv$origPos[1:3] = c('66044','66095','66144')
tsnv$ref[1:3] = c('T','T','T')
tsnv$alt[1:3] = c('C','C','C')

tsnv = attach_gnomad_MAF(tsnv)
tsnv$gnomad_MAF
testthat::test_that("gnomad AF extraction works", {
testthat::expect_equal(gnomad_MAF(chr = '19', position = '66044', ref = 'T', alt = 'C', gnomad = GNOMAD), 0.0245902)
})
