# fastFMM: Fast Functional Mixed Models using Fast Univariate Inference (FUI)

[![](http://cranlogs.r-pkg.org/badges/fastFMM)](https://CRAN.R-project.org/package=fastFMM)

## Repository Description

Repository for the development version of the R Package `fastFMM`. For more information, see the official `fastFMM` $\texttt{CRAN}$ [site](https://CRAN.R-project.org/package=fastFMM).  

## `fastFMM` R Package

## Installation

Download the $\texttt{R}$ Package `fastFMM` by running the following command within $\texttt{R}$ or $\texttt{RStudio}$:

```{R}
install.packages("fastFMM", dependencies = TRUE)
```

Alternatively, the development version of the $\texttt{R}$ Package `fastFMM` can be downloaded as follows:

```{R}
library(devtools)
install_github("gloewing/fastFMM")
```

###  Package Usage

- For usage and a tutorial on package functions, please refer to [fastFMM's Vignette](https://rpubs.com/gloewinger/1110512). 
- For more extensive examples, see the [Photometry Analysis Guide](https://github.com/gloewing/photometry_FLMM/blob/main/README.md)
<br />

## Repository Folders
1) The 'R' folder contains the code of the package, including `fui.R` and `plot_fui.R`. The `plot_fui.R` is still under development and has not been widely tested.

2) The 'vignettes' folder contains a vignette which shows how to use different arguments of the `fui` function. This vignette can also be viewed in the link above (under Package Usage). 

<br />

## Dataset Links

The example data set is available in the 'vignettes' folder under the name 'time_series.csv'.

## Calling fastFMM from Python

See the Python version of [Photometry FLMM Guide Part I](https://github.com/gloewing/photometry_FLMM/blob/main/Tutorials/Photometry%20FLMM%20Guide%20Part%20I/fastFMM-photometry-binary.ipynb) for an example of using `fastFMM` in Python through the Python packages `rpy2` and `fast_fmm_rpy2`. The tutorial assumes the `fastFMM` R package (and all its dependenices), and the `rpy2` Python package have already been installed, but we are working on more documentation for how to install and set up these packages for Python users.  Even if you intend to use the package purely within Python, it may be helpful to first install `fastFMM` from within RStudio to ensure all package dependenices are installed automatically. Finally, see 'python_fastFMM_vignette.py' in the [Github repo](https://github.com/gloewing/photometry_FLMM/tree/main/Tutorials) for a very brief example of using `fastFMM` on Python through the Python package `rpy2`. 
