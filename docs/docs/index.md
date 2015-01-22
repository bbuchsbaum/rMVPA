## rMVPA Documentation 

# Overview

The goal of rMVPA is to provide an easy to use library and set of scripts for conduction MVPA analyses. The primary way in which the tools is intended to be used is through a combination of configuration files and a command line interface. No programming is required to run MVPA analyses with rMVPA. Standard multivariate "searchlight" and region of interest analyses can be run using R scripts that are distributed with the package. rMVPA takes advantage of the [caret](http://topepo.github.io/caret/index.html) package to make a large number of machine learning algorithms available for MVPA analyses. Thus, any machine learning model for which caret provides an interface can be accessed through rMVPA, although we will recommend a set of "tried-and-true" models for general use.


# Installation

Installation is currently done through Github. First you need to install **git** for your operating system. You also need to install a recent version (R > 3.1) of the R programming language. The installation instructions below assume toyu are using a Unix-based operating system (MAc-OSX, Linux, etc.).

Open a terminal. Download the package with git:

```
git clone git@github.com:bbuchsbaum/rMVPA.git
R CMD install rMVPA
```

Next you need to download the command line scripts. You should place these scripts in a folder that is on your PATH. For example, on my slocal computer they are in /home/brad/bin.

```
wget https://raw.githubusercontent.com/bbuchsbaum/rMVPA/master/scripts/MVPA_Searchlight.R
wget https://raw.githubusercontent.com/bbuchsbaum/rMVPA/master/scripts/MVPA_Regional.R
```

To make these files executable:

``` 
chmod +x MVPA_Searchlight.R
chmod +x MVPA_Regional.R
```