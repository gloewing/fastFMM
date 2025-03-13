# fastFMM 0.4.0

* Provided pointers to a Python package to call fastFMM from Python.
* Provided pointers to user guides written in Python.
* Updated reference/citations on documentation.

# fastFMM 0.3.0

* Fixed bugs.
* Added (optional) parallelization of step 3.2 in analytic inference fui(), leading to substantial speed ups of fui().
* Added in parallelization functionality for PCs.
* Added in code to remove rows with missing functional outcome values and added in option to impute with longitudinal FPCA (experimental feature).
* Changed default method of moments estimator to MoM=1 (appears to perform comparably to MoM=2 but is much faster and less memory intensive).
* Removed some fui() arguments that were not in use.
  
# fastFMM 0.2.0

* Fixed several bugs.

# fastFMM 0.1.0

* Initial CRAN submission.
