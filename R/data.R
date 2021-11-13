#' GBIF Occurence Dataset
#'
#' A dataset containing the raw GBIF occurrence information for Central American Amphibians. Obtained via download from GBIF on 11/12/2021.
#' @references  GBIF.org (12 November 2021) GBIF Occurrence Download https://doi.org/10.15468/dl.c8gwtr
#'asdasdasd
#' @format A data frame with 250524 rows and 50 columns
#' \Describe{
#'  \item{ID}{}
#'
#'
#' }
#'
#' @usage data(gbif)
#' @doctype data
#'
#' @examples
#' data(gbif)
#' gbif_clean <- clean_gbif(gbif)
#' paste0(nrow(gbif) - nrow(gbif_clean)," records were cleaned from the raw dataset")
#'
#' @source https://doi.org/10.15468/dl.c8gwtr
#"inst/extdata/gbif.csv.gz"
#'
