## LipidFete

![GitHub release](https://img.shields.io/badge/release-v1.0.0-blue.svg)

*LipidFete* is an R package for identifying tendencies of lipid features

## Description
*LipidFete* is an R package for identifying tendencies of lipid features. *LipidFete* integrates the neighbor information to improve the statistical power and visualize the tendencies.


## Installation
LipidFete v1.0.0 is now available in github
```r
library(devtools)
install_github("Hung-Ching-Chang/LipidFete")
```

## Example
### 2D Lipid feature
```r
library(LipidFete)
data(lipid2D)
X <- t(as.matrix(lipid2D[,-c(1:2)]))
X.info <- lipid2D[,1:2]
group <- rep(c(0, 1), c(52,32))
test.result <- LipidFete.test(X = X,
                              X.info = X.info,
                              group = group,
                              radius = 3,
                              own.contri = 0.5,
                              x.distance = 2,
                              y.distance = 1,
                              dimension = 2,
                              permute.time = 10000)

region.plot.2D(test.result = test.result,
               pval.thres = 0.05,
               x.distance = 2,
               y.distance = 1)
```
Please ignore the warning message, which indicates the number of features with p-values greater than 0.1.

### 1D Lipid feature
```r
library(LipidFete)
data(lipid1D)
X <- t(as.matrix(lipid1D[,2:85]))
X.info <- lipid1D[,1, drop = FALSE]
group <- rep(c(0, 1), c(52,32))
test.result <- LipidFete.test(X = X,
                              X.info = X.info,
                              group = group,
                              radius = 2,
                              own.contri = 0.5,
                              x.distance = 2,
                              y.distance = 1,
                              dimension = 1,
                              permute.time = 10000)
region.plot.1D(X = X,
               test.result = test.result,
               group = group,
               pval.thres = 0.05)

```
## License
This software is licensed under MIT.
