## ----setup, echo = FALSE, eval = FALSE----------------------------------------
#  knitr::opts_chunk$set(comment = "#>", warning=FALSE, message=FALSE)
#  library(fastFMM)
#  output: pdf_document
#  # output: rmarkdown::html_vignette

## ----echo = FALSE-------------------------------------------------------------
# Thanks to Yihui Xie for providing this code
# #  %\VignetteEngine{knitr::rmarkdown} 
# %\VignetteEngine{rmarkdown::render}
library(knitr)
hook_output <- knit_hooks$get("output")
knit_hooks$set(output = function(x, options) {
   lines <- options$output.lines
   if (is.null(lines)) {
     return(hook_output(x, options))  # pass to default hook
   }
   x <- unlist(strsplit(x, "\n"))
   more <- "..."
   if (length(lines)==1) {        # first n lines
     if (length(x) > lines) {
       # truncate the output, but add ....
       x <- c(head(x, lines), more)
     }
   } else {
     x <- c(more, x[lines], more)
   }
   # paste these lines together
   x <- paste(c(x, ""), collapse = "\n")
   hook_output(x, options)
 })

## ---- eval=FALSE--------------------------------------------------------------
#  library(devtools)
#  install_github("gloewing/fastFMM")

## ---- eval = FALSE------------------------------------------------------------
#  # install packages (will only install them for the first time)
#  list.of.packages = c("lme4", "parallel", "cAIC4", "magrittr","dplyr",
#                        "mgcv", "MASS", "lsei", "refund","stringr", "Matrix", "mvtnorm",
#                        "arrangements", "progress", "ggplot2", "gridExtra", "here")
#  new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#  if(length(new.packages)){
#    chooseCRANmirror(ind=75)
#    install.packages(new.packages, dependencies = TRUE)
#  }else{
#    # load packages if already installed
#    library(lme4)
#    library(parallel)
#    library(cAIC4)
#    library(magrittr)
#    library(dplyr)
#    library(mgcv)
#    library(MASS)
#    library(lsei)
#    library(refund)
#    library(stringr)
#    library(Matrix)
#    library(mvtnorm)
#    library(arrangements)
#    library(progress)
#    library(ggplot2)
#    library(gridExtra)
#  }

## -----------------------------------------------------------------------------
library(fastFMM) # load our package

## -----------------------------------------------------------------------------
dat <- read.csv("time_series.csv") # read in data
head(dat[,1:6])

## ---- eval = FALSE------------------------------------------------------------
#  Y_mat <- dat[,-seq(1,3)]
#  head(Y_mat[,1:5])

## ---- eval = FALSE------------------------------------------------------------
#  dat <- data.frame(Y = Y_mat, dat[,seq(1,3)])

## ---- eval = FALSE------------------------------------------------------------
#  mod <- fui(Y ~ treatment + # main effect of cue
#                    (1 | id),  # random intercept
#                  data = dat)

## ---- eval = FALSE------------------------------------------------------------
#  mod <- fui(Y ~ treatment + # main effect of cue
#                    (treatment | id),  # random intercept
#                    data = dat,
#                    parallel = TRUE,
#                    analytic = FALSE) # bootstrap

## ---- eval = FALSE------------------------------------------------------------
#  Y_mat <- dat[,-seq(1,3)])
#  L <- ncol(Y_mat) # number of columns of functional outcome
#  
#  mod <- fui(Y ~ treatment + # main effect of cue
#                    (treatment | id),  # random intercept
#                    data = dat,
#                    argvals = seq(from = 1, to = L, by = 3) # every 3rd data point
#                    )

## ---- eval = FALSE------------------------------------------------------------
#  Y_mat <- dat[,-seq(1,3)])
#  L <- ncol(Y_mat) # number of columns of functional outcome
#  
#  # model 1: random slope model
#  mod1 <- fui(Y ~ treatment + # main effect of cue
#                    (treatment | id),  # random intercept
#                    data = dat,
#                    var = FALSE)
#  
#  # model 2: random intercept model
#  mod2 <- fui(Y ~ treatment + # main effect of cue
#                    (1 | id),  # random intercept
#                    data = dat,
#                    var = FALSE)
#  
#  # compare model fits
#  colMeans(mod1$aic)
#  colMeans(mod2$aic)

## ---- eval = FALSE------------------------------------------------------------
#  mod <- fui(Y ~ treatment + # main effect of cue
#                    (treatment | id/trial),  # random intercept
#                    data = dat,
#                    subj_ID = "id")

## -----------------------------------------------------------------------------
mod <- fui(Y ~ treatment + # main effect of cue
                  (treatment | id),  # random intercept
                  data = dat) 

fui_plot <- plot_fui(mod)

## -----------------------------------------------------------------------------
align_time <- 1 # cue onset is at 2 seconds
sampling_Hz <- 15 # sampling rate
# plot titles: interpretation of beta coefficients
plot_names <- c("Intercept", "Mean Signal Difference: Cue1 - Cue0") 
fui_plot <- plot_fui(mod, # model fit object
                     x_rescale = sampling_Hz, # rescale x-axis to sampling rate
                     align_x = align_time, # align to cue onset
                     title_names = plot_names,
                     xlab = "Time (s)",
                     num_row = 2) 

## -----------------------------------------------------------------------------
head(fui_plot$treatment)

