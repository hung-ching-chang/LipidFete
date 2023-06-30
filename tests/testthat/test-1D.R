# Load example data
data(lipid1D)
X <- t(as.matrix(lipid1D[,2:85]))
X.info <- lipid1D[,1]
group <- rep(c(0, 1), c(52,32))

# start tests
test_that("region.plot.1D() function computed successfully.",{
  expect_type(X.info, "integer")
  expect_type(group, "double")
  test.result <- LipidFete.test(X = X,
                                X.info = X.info,
                                group = group,
                                radius = 2,
                                own.contri = 0.5,
                                x.distance = 2,
                                y.distance = 1,
                                dimension = 1,
                                permute.time = 10000)
  expect_type(test.result$direction, "character")
  expect_type(test.result$smoothing.pval.BH, "double")
  expect_type(test.result$marginal.pval.BH, "double")
  expect_type(test.result$log2.FC, "double")
  figure.1D <- region.plot.1D(X = X,
                              X.info = X.info,
                              group = group,
                              direction = test.result$direction,
                              smoothing.pval = test.result$smoothing.pval.BH,
                              marginal.pval = test.result$marginal.pval.BH,
                              feature.name = colnames(X)[1],
                              cut.point = 0.05)
  expect_s3_class(figure.1D, "ggplot")
})
