# PatternMatchR: A Tool for Matching and Visualiszation of Omics Patterns

PatternMatchR is a tool for finding and visualization of matching patterns in Omics.

The PatternMatchR package is designed to find and visualize matching patterns in Omics analyses. It supports the user by providing
- Manhattan and volcano plots for CpG selection;
- Trait Methylation plots;
- Methylation profile plots and
- Correlation plots.

## Installation

library('remotes')  
install_github('SteRoe/PatternMatchR')

## Usage

Point your working directory to the location with your files to analyze. This folder should also contain a valid yml-file with all necessary settings for

dataDir1: a list containg Omics results for a first matching scenario

dataDir2: a list containg Omics results for a second matching scenario

dataDir3: a list containg Omics results for a third matching scenario

workDir: working directory

P_VALWarningThreshold: unusual low p-value for which a warning should be generated (5e-50)

debugMode: TRUE


An example is given in the ./inst/extdata/ folder.
Detailed explanations of data structures are in the vignette.

setwd("your.working.directory")  
Start the App using: PatternMatchR::PatternMatchRApp()

## License
[MIT](https://choosealicense.com/licenses/mit/)
