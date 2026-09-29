# install.R — R dependencies for the forecast skill score analysis
# Run once:  Rscript install.R

cran <- c(
  "arrow", "tidyverse", "forecast", "boot", "zoo",
  "MMWRweek", "fitdistrplus", "foreach", "doParallel",
  "epidatr", "httr", "xml2", "remotes"
)
install.packages(cran, repos = "https://cloud.r-project.org")

# Not on CRAN — reichlab GitHub packages (used only by the data-collection step)
remotes::install_github("reichlab/covidData")
remotes::install_github("reichlab/covidHubUtils")
