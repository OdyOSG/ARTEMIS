#Cannot test getCohortSet due to its reliance on an active cdm connection

test_that("stringDF_from_cdm works", {

  fakeConDF <- as.data.frame(matrix(c("1","1","2","3",
                                      "2021-09-06","2021-12-06","2022-09-06","2020-09-06",
                                      "1338512","1338512","1338512","1",
                                      "1338512","1338512","1338512","1",
                                      "doxorubicin","doxorubicin","doxorubicin","INVALID_DRUG"),
                                    ncol = 5))

  colnames(fakeConDF) <-c("person_id","drug_exposure_start_date","drug_concept_id",
                          "ancestor_concept_id","concept_name")

  fakeDrugs <- as.data.frame(matrix(c("Doxorubicin","35803022",
                                      "35803022","Doxorubicin",
                                      "1338512","Drug","Ingredient",
                                      ""), ncol = 8, nrow = 1))

  colnames(fakeDrugs) <- c("Raw_HemOnc","concept_id","Manual","concept_me",
                           "concept_id_2","domain_id","concept_class_id",
                           "Manual_Req")

  stringDF <- stringDF_from_cdm(con_df = fakeConDF, writeOut = F, outputName = "Output", validDrugs = fakeDrugs)

  expect_identical(stringDF$person_id,c("1","2"))
  expect_identical(stringDF$seq[1],"0.doxorubicin;91.doxorubicin;")
  expect_identical(stringDF$seq[2],"0.doxorubicin;")

})

test_that("filter_stringDF works", {
  fakeStringDF <- as.data.frame(matrix(c("Test1","Test2","0.A;1.B","0.A;1.B;0.A;1.B",2,4),ncol=3))

  colnames(fakeStringDF) <- c("person_id","seq","Count")

  filteredStringDF <- fakeStringDF %>% filter_stringDF(min = 3)

  expect_equal(dim(filteredStringDF)[1],1)

})

