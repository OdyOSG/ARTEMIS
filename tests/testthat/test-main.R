test_that("generateRawAlignment works", {
  fakeStringDF <- as.data.frame(matrix(c("Test1","Test2",
                                         "0.A;7.C;0.A;7.C","0.A;7.B;0.A;7.B",
                                         4,4),ncol=3))

  colnames(fakeStringDF) <- c("person_id","seq","Count")

  fakeRegimens <- as.data.frame(matrix(c(1,2,
                                          "Test","Test",
                                          "Test","Test",
                                          1,1,
                                          "Test","Test",
                                          1,1,
                                          "Reg1","Reg2",
                                          1,1,
                                          1,2,
                                          "AB","AC",
                                          "1,7","1,7",
                                          "","",
                                          7,7,
                                          1,1,
                                          "","",
                                          "","",
                                          "TRUE","TRUE",
                                          "Cont.","Cont.",
                                          "0.A;7.C;","0.A;7.B;",
                                          "0.a;7.c;","0.a;7.b;"),ncol=20))

  colnames(fakeRegimens) <- c("regCodeExt","metaCondition","condition",
                              "conditionCode","context","contextCode",
                              "regName","variant","regCode","component",
                              "day","cycleTaken","cycleLength","noCycles",
                              "branchInfo","Radio.Therapy.","continuous",
                              "noCycles_Original","regString","shortString")

  fakeRawAlignments <- generateRawAlignments(stringDF = fakeStringDF,
                                             regimens = fakeRegimens,
                                             g = 0.4,
                                             Tfac = 0.5,
                                             verbose = 0,
                                             mem = -1,
                                             removeOverlap = 1,
                                             method = "PropDiff",
                                             writeOut = F, outputName = "test")

  unlink(here::here("output/test.csv"))

  #Scores
  expect_equal(fakeRawAlignments$Score, c("","1.6","1.5625",
                                          "","1.6","1.5625"))
  #Regimens
  expect_equal(fakeRawAlignments$regName, c("Reg1","Reg1","Reg1",
                                            "Reg2","Reg2","Reg2"))

  #Patients
  expect_equal(fakeRawAlignments$personID, c("Test1","Test1","Test1",
                                             "Test2","Test2","Test2"))

})

test_that("processAlignments works", {
  fakeResults <- as.data.frame(matrix(c("Reg1","Reg1","Reg1",
                                        "0.A;7.A","0.A;7.A","0.A;7.A",
                                        "0.A;7.A;7.A;","0.A;7.A","7.A;7.A",
                                        "","2","2",
                                        "","1","1",
                                        "","2","2",
                                        "","1","2",
                                        "","2","3",
                                        "","2","2",
                                        "","2","2",
                                        NA,1,1,
                                        "Test1","Test1","Test1"), ncol = 12))

  colnames(fakeResults) <- c("regName","Regimen","DrugRecord","Score","regimen_Start",
                             "regimen_End","drugRec_Start","drugRec_End",
                             "Aligned_Seq_len","totAlign","adjustedS","personID")

  processedResults <- processAlignments(fakeResults, writeOut = F,
                                        regimenCombine = 28, outputName = "test")

  unlink(here::here("output/test.csv"))

  expect_equal(processedResults$t_start, 0)
  expect_equal(processedResults$t_end, 14)
  expect_equal(processedResults$adjustedS, 1)

})
