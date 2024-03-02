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

For the usage and a tutorial on package functions, please refer to [fastFMM's Vignette](https://rpubs.com/gloewinger/1110512). 

<br />

## Repository Folders
1) The 'R' folder contains the code of the package, including `fui.R` and `plot_fui.R`. The `plot_fui.R` is still under development and has not been widely tested.

2) The 'vignettes' folder contains a vignette which shows how to use different arguments of the `fui` function. This vignette can also be viewed in the link above (under Package Usage). 

<br />

## Dataset Links

The example data set is available in the 'vignettes' folder under the name 'time_series.csv'.

## Calling fastFMM from Python

See 'python_fastFMM_vignette.py' in the [Github repo](https://github.com/gloewing/photometry_FLMM/tree/main/Tutorials) for a brief example of using `fastFMM` on Python through the Python package `rpy2`. We are working on more documentation. The tutorial assumes the `fastFMM` R package (and all its dependenices), and the `rpy2` Python package have already been installed. Even if you intend to use the package purely within Python, it may be helpful to first install `fastFMM` in RStudio to ensure all package dependenices are installed automatically.
