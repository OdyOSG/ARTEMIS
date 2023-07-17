test_that("align works", {

  g <- 0.4
  Tfac <- 0.25
  verbose = 0
  mem = as.integer(-1)
  removeOverlap = 1

  regimen1 <- encode("0.carboplatin;14.cisplatin;4.Doxorubicin")
  drugRecord <- encode("0.Carboplatin;10.Leuprolide;4.Cisplatin;4.Doxorubicin")

  output <- align(regimen1,"Test",drugRecord,g,Tfac,NA,verbose,mem,removeOverlap,"PropDiff")

  expect_equal(output$adjustedS[2], 0.86666667)
})

test_that("align errors1", {

  g <- 0.4
  Tfac <- 0.25
  verbose = 0
  mem = as.integer(-1)
  removeOverlap = 1

  regimen1 <- encode("0.carboplatin;14.cisplatin;4.Doxorubicin")
  regimen2 <- encode("0.carboplatin;7.cisplatin;7.Doxorubicin")
  drugRecord <- encode("0.Carboplatin;10.Leuprolide;4.Cisplatin;4.Doxorubicin")

  expect_output(align(list(regimen1,regimen2),"Test",drugRecord,g,Tfac,NA,verbose,mem,removeOverlap,"PropDiff"),
                "Multiple regimens but only one regname. Please check regnames.")
})

test_that("align works with multiple inputs", {

  g <- 0.4
  Tfac <- 0.25
  verbose = 0
  mem = as.integer(-1)
  removeOverlap = 1

  regimen1 <- encode("0.carboplatin;14.cisplatin;4.Doxorubicin")
  regimen2 <- encode("0.carboplatin;7.cisplatin;7.Doxorubicin")
  drugRecord <- encode("0.Carboplatin;10.Leuprolide;4.Cisplatin;4.Doxorubicin")

  output <- align(list(regimen1,regimen2),list("Reg1","Reg2"),drugRecord,g,Tfac,NA,verbose,mem,removeOverlap,"PropDiff")

  expect_equal(output$regName,c("Reg1","Reg1","Reg2","Reg2"))
  expect_equal(output$adjustedS,c(NA,0.8666666667,NA,0.7892857143))

})

