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
