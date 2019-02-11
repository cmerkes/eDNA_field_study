# Code to assist the USFWS with eDNA field sampling designs for eDNA


Richard A. Erickson, Christopher M. Merkes, and Erica L. Mize


This code supports the manuscript "Refinement of eDNA as an early monitoring tool at the landscape-level: Study design considerations".
The code fit a three-level occupancy model with two detection assays in the program [Stan](mc-stan.org) through [R](https://www.r-project.org/).
We programmed this code for a specific project with the USFWS and it will likely not generalize to other situations.
The primary purpose of releasing this code is transparency and scientific reproducible. 
However, this code could be adapted by people who understand Stan and want to create their own three-level occupancy model.

## Code files

This repository contains the following files and folder:
- `README.md`: This file
- `LICENSE`: The standard USGS software license
- `model_with_coef.R`: The R Script that formats data and runs the Stan model
- `coefModel`: A folder with script
  - `positiveSampleCoef.stan`: The Stan model that is called by the R script.
  - `simSampleSize.R`: R script that takes the output from the model and calculates the number of samples required for one positive detection.

## Data 

The data used in this project belong to the USFWS and the point of contact is Erica Mize: erica_mize@fws.gov


## Contact for code 

Primary code developer:  Richard A. Erickson (rerickson@usgs.gov)

## Disclaimer

This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey, an agency of the United States Department of Interior. For more information, see the [official USGS copyright policy](https://www2.usgs.gov/visual-id/credit_usgs.html#copyright/).


This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.

This software is provided "AS IS".
