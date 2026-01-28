## Test environments

* local Ubuntu 22.04.3, R 4.1.2
* Fedora Linux, R-devel, clang, gfortran (on R-hub)
* Ubuntu Linux 20.04.1 LTS, R-release, GCC (on R-hub)
* Windows Server 2022, R-devel, 64 bit (on R-hub)
* Mac OS 13.3, R-release, clang (on mac.r-project.org)


## R CMD check results

1. The following note showed up consistently on several testing platforms:
```
* checking installed package size ... NOTE
  installed size is 13.4Mb
  sub-directories of 1Mb or more:
    libs  12.9Mb
```
The exact size of the installed package and libraries varied, depending on
the platform. It is my understanding that such size is not unexpected in
the case of packages with compiled code

I think this package was archived because a re-submission  was submitted with the
same release (0.6.0) as the release that I was trying to correct (0.6.0). This 
submission (0.6.1) attempts to correct this in line with CRAN policy (https://cran.r-project.org/web/packages/policies.html).


2. Two that are only found on Windows (Server 2022, R-devel 64-bit): 

```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```
As noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX and can likely be ignored.

3. The third is 

```
* checking for non-standard things in the check directory ... NOTE
Found the following files/directories:
  ''NULL''
```

As noted in [R-hub issue #560](https://github.com/r-hub/rhub/issues/560), this seems to be an Rhub issue and so can likely be ignored. 

4. A fourth that is found with *Fedora Linux, R-devel, clang, gfortran* and *Ubuntu Linux 20.04.1 LTS, R-release, GCC*

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
```

This also seems to be a recurring issue on Rhub [R-hub issue #560](https://github.com/r-hub/rhub/issues/548) and so can likely be ignored.

---

This version restores dependency on the Boost library. The previous version was using the implementation of certain functions directly from the std library instead, but these functions are not currently supported by the clang compiler.

## Downstream dependencies
GeoModels, VeccTMVN