### MP Variant Table
source("/home/ionadmin/github_app/GDrive_VariantReport/Gauths.R")
# mpvarlink = googledrive::as_id("https://docs.google.com/spreadsheets/d/1B-NfpRNhadl7w1f5UPkRA_XEg4YI3N4pHRxd9yZgZkc/edit?usp=drive_web&ouid=116704210424700639172")
MPvars = googlesheets4::read_sheet(mpvarlink, skip = 1)
mpvs = MP_check_retrieve_table(snvt, MPvars)
mpv_filepath = paste0(analysis_dir,'/annotation_output/annotation_MP_variant.tsv')
if(nrow(mpvs) > 0){
  readr::write_tsv(mpvs, file = mpv_filepath )
}


### BRAF Variant Class
# braflink = googledrive::as_id("https://docs.google.com/spreadsheets/d/1xQ3FfHV2JLndT7J_yLHGcrb7Uqpc-cjjio9OYXkjz14/edit")
BRAFinfo = googlesheets4::read_sheet(braflink, skip = 1)
brafi = BRAF_class_retrieve_table(snvt, BRAFinfo)
braf_class = paste0(analysis_dir,'/annotation_output/annotation_BRAF_variant_class.tsv')
if(nrow(brafi) > 0){
  readr::write_tsv(brafi, file = braf_class)
}

### HRR genes
# hrrgenes = googledrive::as_id('https://docs.google.com/spreadsheets/d/11FZz5m34IYmK1-o2UkjFJ7gatXUBr4_F0mhVdGMxQm8/edit#gid=230169235')
HRRgenes = googlesheets4::read_sheet(ss = hrrgenes, sheet  = "HRR_single_genes", skip = 1)
hrrs = HRR_check_retrieve_table(snvt, HRRgenes)
HRRgenepath = paste0(analysis_dir,'/annotation_output/annotation_HRRgenes.tsv')
if(nrow(mpvs) > 0){
  readr::write_tsv(hrrs, file = HRRgenepath)
}


## External filter chains
## https://docs.google.com/spreadsheets/d/16qltcg5StmQk0A2Spb6iBOSXkSgZPqQ_3onqAFQt5tI/edit#gid=0

