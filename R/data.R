#' gbif
#'
#' GBIF Occurence Dataset
#'
#' A dataset containing the raw GBIF occurrence information for Central American Amphibians. Obtained via download from GBIF on 11/12/2021.
#' @references  GBIF.org (12 November 2021) GBIF Occurrence Download https://doi.org/10.15468/dl.c8gwtr
#'
#' @format A data frame with 250524 rows and 50 columns
#' \Describe{
#'  \item{ID}{}
#'
#'
#' }
#'
#' @usage gbif <- read.csv(system.file('extdata','gbif.csv.gz', package = "evogear"))
#'
#' @examples
#' data(gbif)
#' gbif_clean <- clean_gbif(gbif)
#' paste0(nrow(gbif) - nrow(gbif_clean)," records were cleaned from the raw dataset")
#'
#' @source https://doi.org/10.15468/dl.c8gwtr


#'iucn
#'
#'IUCN Occurence dataset
#'
#'A dateset containing amphibian species in Central America. Downloaded from IUCN.
#'Filtered based on species found in filtered GBIF dataset (IUCN <- IUCN_maps[which(IUCN_maps$binomial %in% pp$species),])
#'
#'@format A .shp file that should be loaded as a spdf with 437 elements with 27 variables
#'@usage ReadOGR('inst/extdata/iucn.shp')




