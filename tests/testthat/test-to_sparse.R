test_that("to_sparse works", {

  obj <- data.frame(gene=paste0("gene",seq_len(10)),
                    matrix(data = 0,nrow=10,ncol = 20))
  obj2 <- to_sparse(obj)
  obj3 <- to_sparse(obj = obj2)
  obj4 <- to_sparse(as.matrix(obj2))
  testthat::expect_true(methods::is(obj2,"sparseMatrix"))
  testthat::expect_true(methods::is(obj3,"sparseMatrix"))
  testthat::expect_true(methods::is(obj4,"sparseMatrix"))
})
