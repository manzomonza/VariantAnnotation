# ## 0. Adjust protein three letter to one letter
amino_acid_code_three_to_one = function(snv){
  snv$protein = unname(sapply(snv$protein, amino_acid_conversion_three_to_one))
  return(snv)
}
