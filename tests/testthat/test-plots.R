test_that("plotProcessed works", {

  fakePlot <- as.data.frame(matrix(c(10,20,35,
                                     19,29,49,
                                     "Test1","Test2","Test2",
                                     "Yes","Yes","Yes",
                                     1,1,0.75),
                                   ncol = 5))

  colnames(fakePlot) <- c("t_start","t_end","component","regimen","adjustedS")

  fakePlot$t_start <- as.numeric(fakePlot$t_start)
  fakePlot$t_end <- as.numeric(fakePlot$t_end)

  plotProcesssed(fakePlot, fontSize = 2.5)

  suppressWarnings(

    vdiffr::expect_doppelganger("testPlot",plotProcesssed(fakePlot, fontSize = 2.5))

  )

})

test_that("plotOutput works with complex output", {

  fakeResults <- as.data.frame(matrix(c("Reg1","Reg1","Reg1","Reg1","Reg2","Reg3","Reg4","Reg5",
                                        "0.A;7.A","0.A;7.A","0.A;7.A","0.A;7.A","0.A;3.A;3.A","0.A;3.A","1A.1A","1A.1B",
                                        "0.A;7.A;7.A;28.A;7.A;1.A;1.B","0.A;7.A","7.A;7.A","28.A;7.A","7.A;7.A","7.A;7.A","1.A;1.A","1.A;1.A",
                                        "","2","2","2","1.2","1","1","1",
                                        "","1","1","1","1","1","2","2",
                                        "","2","2","2","3","2","2","2",
                                        "","0","2","4","2","4","6","6",
                                        "","2","3","6","3","5","7","7",
                                        "","2","2","2","2","3","2","2",
                                        "","2","2","2","2","2","2","2",
                                        NA,1,1,1,0.6,0.5,0.75,0.75,
                                        "Test1","Test1","Test1","Test1","Test1","Test1","Test1","Test1"), ncol = 12))

  colnames(fakeResults) <- c("regName","Regimen","DrugRecord","Score","regimen_Start",
                             "regimen_End","drugRec_Start","drugRec_End",
                             "Aligned_Seq_len","totAlign","adjustedS","personID")

  output <- plotOutput(output = fakeResults, fontSize = 2.5,
                       regimenCombine = 1, returnDat = T)

  expect_equal(output$t_start,c(0,42))
  expect_equal(output$t_end,c(14,50))
  expect_equal(output$component,c("Reg1","Reg1"))

})
