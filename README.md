# GPvecchia
Fast Gaussian-process inference using general Vecchia approximations

For examples of how to use the package, please see the vignettes folder. Please note that GPvecchia is under active development and not stable at this time.

## References
[Katzfuss, M., & Guinness, J. (2017). A general framework for Vecchia approximations of Gaussian processes. *arXiv:1708.06302*.](https://arxiv.org/abs/1708.06302)

[Katzfuss, M., Guinness, J., & Gong, W. (2018). Vecchia approximations of Gaussian-process predictions. *arXiv:1805.03309*.](https://arxiv.org/abs/1805.03309)

## Installation
To ensure that the algorithms run efficiently on your computer, we recommend installing the package by downloading the repo to a local folder, setting your R working directory to that folder, and then running the following code in R:
```{r}
library(GpGp); library(Matrix); library(RcppParallel)
library(parallel); library(sparseinv); library(fields)
for (nm in list.files('GPvecchia/R',pattern = "\\.[RrSsQq]$")) {
  cat(nm,":"); source(file.path('GPvecchia/R',nm)); cat("\n")
}
Rcpp::sourceCpp('GPvecchia/src/U_NZentries.cpp')
Rcpp::sourceCpp('GPvecchia/src/MaxMin.cpp')
```

<!--- 
This package can be installed directly from R by running the following code:
```{r}
library(devtools)
install_github("katzfuss-group/GPvecchia")
```
Alternatively, one can download the .tar.gz file from the main directory here and then run:
```{r}
install.packages("GPvecchia_0.1.tar.gz", repos = NULL, type = "source")
```
-->

Note that [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is required for compiling C/C++ with OpenMP on Windows systems. When installing Rtools, the system PATH needs to be set so that the C++ compiler included in Rtools can be found by R. Once Rtools is installed, `system('g++ -v')` can be used to check if the compiler is accessible from R.

