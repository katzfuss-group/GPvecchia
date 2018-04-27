# GPvecchia
Fast Gaussian-process inference using general Vecchia approximations

## Installation
This package can be installed directly from R by running the following code:
```{r}
library(devtools)
install_github("katzfuss-group/GPvecchia")
```
Alternatively, one can download the .tar.gz file from the main directory here and then run:
```{r}
install.packages("~/GPvecchia_0.1.tar.gz", repos = NULL, type = "source")
```

Note that [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is required for compiling C/C++ with OpenMP on Windows system. When installing Rtools, the system PATH needs to be set so that the C++ compiler included in Rtools can be found by R. Once Rtools is installed, `system('g++ -v')` can be used to check if the compiler is accessible from R.

