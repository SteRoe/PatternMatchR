# PatternMatchR: A Tool for Matching and Visualiszation of Omics Patterns

PatternMatchR is a tool for finding and visualizing of matching patterns in Omics.

The PatternMatchR package is designed to find and visualize matching patterns in Omics analyses. Examples are for epigenetics. 

It supports the user by providing
- select and merge data out of various data sources;
- reduce number of results by defining limits for p-values;
- reduce number of traits by clustering similar traits together;
- plotting a heatmap which shows similarities between traits. Data from up to three different topics are marked with red/green/blue colors inside the heatmap in order to easily visualize data source of a certain trait;
- plotting a scatter plot matrix (SPLOM), which shows the relationship between similar traits and methylation.

## Installation

library('remotes')  
install_github('SteRoe/PatternMatchR')

## Usage
PatternMatchR is configured using a <config.yml> file inside the main folder (working directory of the application). The structure of this configuration file is described in the application vignette.

An example for used data files is given in the ./inst/extdata/ folder.
Detailed explanations of data structures are given in the vignette.

setwd("your.working.directory")  
Start the App using: PatternMatchR::PatternMatchRApp()

## License
[MIT](https://choosealicense.com/licenses/mit/)
