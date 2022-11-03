#' Collapse shapefiles in an sf object based on unique values in a column
#'
#'@param sf the sf object you wish to collapse
#'@param group_col the grouping column
#'
#'
#'@export
#'

sf_collapse <- function(sf, group_col){

  sf[is.na(sf)] <- "NA"

  func <- function(x){
    x%>%
      st_combine()%>%
      st_union()
  }

  #split by category and union it
  sf %>%
    split(., f = .[[group_col]]) %>%
    lapply(st_union) -> out

  #grab labels here
  labels <- names(out)

  do.call(c, out) -> out

  #create sf dataframe
  out <- st_as_sf(data.frame(labels,st_geometry(out)))

  #Add name as group_col input
  names(out) <- c(paste(group_col), "geometry")

  return(out)}