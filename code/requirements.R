## This file install necessary packages for scripts accompanying 
## [Introducing Sample Recycling Method]
################################################################################

# packages of interest for this project
packages <- c("parallel", # for parallel processing
              "ggplot2", # for plotting
              "latex2exp") # for LaTeX expressions in plots
              # "knitr", # for basic features to generate LaTeX tables
              # "kableExtra") # for advanced features to generate LaTeX tables 
# extract packages installed 
packages.installed <- as.character(installed.packages()[,1])
# install needed packages 
for (package in packages) {
  if (!(package %in% packages.installed)) {
    install.packages(package, dependencies = TRUE)
  }
}
# remove introduced variables 
rm(package, packages, packages.installed)
