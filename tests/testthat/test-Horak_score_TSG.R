exampletsg = "./testing/annotation_TSG.tsv"

#
resulttable = Horak_score_TSG(exampletsg)
expected_scores = c(0,0,0,8,0)
test_that("Hscore retrieval for TSG or non TSG genes works accurately", {
for(i in seq_along(resulttable$TSG_hscore)){
  expect_equal(resulttable$TSG_hscore[i], expected_scores[i])
  }
})
