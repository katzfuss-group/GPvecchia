# GPvecchia
Fast Gaussian-process inference using Vecchia approximations

## Installation
This package can be installed directly from R with the following
```{r}
library(devtools)
install_github("katzfuss-group/GPvecchia")
```
Note that [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is required for compiling C/C++ with OpenMP on Windows system. Attention needs to be paid in a step where you can edit the system PATH so that the C++ compiler that is included in Rtools can be found by R. Once Rtools is installed, you can use `system('g++ -v')` to check if the compiler is accesible from R.

You can also mannually intall the package in R with .tar.gz file using
`install.packages("~/GPvecchia_0.1.tar.gz", repos = NULL, type = "source")`
