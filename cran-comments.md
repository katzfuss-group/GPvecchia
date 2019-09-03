## Test environments
* local ubuntu 18.04, R 3.6.1
* ubuntu 16.04 (Travis CI), R 3.6.1
* local OS X install, R 3.5.1

## R CMD check results
There were no ERRORs or WARNINGs under either Ubuntu version

There were 2 NOTEs on Travis:

* checking installed package size ... NOTE
  installed size is  5.1Mb
  sub-directories of 1Mb or more:
    libs   4.6Mb

This is not unexpected in packages with compiled code as observed by Dirk Eddelbuettel
https://stackoverflow.com/questions/53819970/r-package-libs-directory-too-large-after-compilation-to-submit-on-cran

* checking for hidden files and directories ... NOTE
  Found the following hidden files and directories:
  	.travis.yml
  These were most likely included in error. See section ‘Package
  structure’ in the ‘Writing R Extensions’ manual.

The said file is listed in .Rbuildignore and the note goes away in the local check of the build tarball


There was 1 NOTE on local ubuntu installation:

* checking installed package size ... NOTE
  installed size is  6.3Mb
  sub-directories of 1Mb or more:
    libs   5.9Mb

See explanation above



There was 1 WARNING under local OS X:

*  checking whether package ‘GPvecchia’ can be installed (41.2s)
   Found the following significant warnings:
     Warning: package ‘GpGp’ was built under R version 3.5.2

Local R version used for development is 3.5.1



## Downstream dependencies
none