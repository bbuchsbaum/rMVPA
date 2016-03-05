
rMVPA
=====

Multivoxel Pattern Analysis in R

### This package is under development.

### Introduction

`rMVPA` is an R library for multivariate pattern analysis of neuroimaging data. The goal of this library is to make MVPA analyses easy. It can be used both programmatically from within R or using a command line interface. 'rMVPA' leverages the 'caret' library for the underlying machine learning interface. What `rMVPA` provides is the infrastructure for conducting machine learning analyses on neuroimaging data. 

### Installation

## Using devtools

To install `rMVPA` from within R, use the `devtools` function `install_github`.

From within R:

```
library(devtools)
install_github("bbuchsbaum/rMVPA")
```

## Using `git` from the command line

```
git clone git@github.com:bbuchsbaum/rMVPA.git
R CMD install rMVPA
```
