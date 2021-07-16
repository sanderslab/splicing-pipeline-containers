install.packages("devtools", repos='http://cran.us.r-project.org')
devtools::install_github("davidaknowles/leafcutter/leafcutter")

library(leafcutter)

args <- commandArgs(trailingOnly=TRUE) 
print(args) 
