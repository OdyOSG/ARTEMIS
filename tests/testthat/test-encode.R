test_that("encode works", {
  expect_equal(encode("1.A;1.B"), list(c("1","A"),c("1","B")))
})
