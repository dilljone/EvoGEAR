#' creates a simple map based on supplied country name list
#'
#' @param country_names A vector of country names to use
#' @export
mapcreate <- function(country_names){

  require(maps)
  require(maptools)

  map <- maps::map('world',
                   regions = country_names,
                   fill = TRUE,
                   interior = FALSE)


  IDs <- sapply(strsplit(map$names,":"),function(x) x[1])

  map_sp <- maptools::map2SpatialPolygons(map, IDs = IDs)%>%
    raster::aggregate(., dissolv = TRUE)

  return(map_sp)

}
