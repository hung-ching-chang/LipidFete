# Load example data
data(lipid2D)
X <- t(as.matrix(lipid2D[,-c(1:2)]))
X.info <- lipid2D[,1:2]
group <- rep(c(0, 1), c(52,32))

# start tests
test_that("LipidFete.test() function computed successfully.",{
  expect_s3_class(X.info, "data.frame")
  expect_type(group, "double")
  test.result <- LipidFete.test(X = X,
                                X.info = X.info,
                                group = group,
                                radius = 3,
                                own.contri = 0.5,
                                x.distance = 2,
                                y.distance = 1,
                                dimension = 2,
                                permute.time = 10000)
  expect_s3_class(test.result, "data.frame")
  expect_equal(nrow(test.result), ncol(X))
  expect_error(LipidFete.test(X = X,
                              X.info = X.info,
                              group = group,
                              radius = 3,
                              own.contri = 1,
                              x.distance = 2,
                              y.distance = 1,
                              dimension = 2,
                              permute.time = 10000))
})


test_that("region.plot.2D() function computed successfully.",{
  test.result <- LipidFete.test(X = X,
                                X.info = X.info,
                                group = group,
                                radius = 3,
                                own.contri = 0.5,
                                x.distance = 2,
                                y.distance = 1,
                                dimension = 2,
                                permute.time = 10000)
  expect_type(test.result$direction, "character")
  expect_type(test.result$smoothing.pval.BH, "double")
  expect_type(test.result$marginal.pval.BH, "double")
  expect_type(test.result$log2.FC, "double")
  figure.2D <- region.plot.2D(X.info = X.info,
                              direction = test.result$direction,
                              smoothing.pval = test.result$smoothing.pval.BH,
                              marginal.pval = test.result$marginal.pval.BH,
                              log2.FC = test.result$log2.FC,
                              cut.point = 0.05,
                              x.distance = 2,
                              y.distance = 1)
  expect_s3_class(figure.2D, "ggplot")
})

test_that("pval.annotation() function computed successfully.",{
  input <- NA
  expect_type(NA, "logical")
  output <- pval.annotation(input)
  expect_type(output, "logical")
})
