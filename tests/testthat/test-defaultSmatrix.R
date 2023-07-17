test_that("defaultSmatrix works", {
  s1 <- encode("1.A;1.B")
  s2 <- encode("1.A;1.C")

  default <- defaultSmatrix(s1,s2)

  testMat <- matrix(c(1,-1.1,-1.1,-1.1,1,-1.1,-1.1,-1.1,1), ncol = 3)

  colnames(testMat) <- c("A","B","C")
  rownames(testMat) <- c("A","B","C")

  expect_equal(as.matrix(default), as.matrix(testMat))
})
