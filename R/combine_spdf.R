#' Combine a list of spdfs together and group by binomialk name
#'
#'@param input a list of lists containing each spdf
#'@param binomial Variable relating to the name of the column where binomials
#'
#'@examples
#'    x <- pp
#' y <- IUCN
#' y <- y[y@data$legend == 'Extant (resident)',]
#' y <- raster::aggregate(y, by = 'binomial',dissolve = TRUE)
#' x$binomial <-  x$species
#'
#'list <- list()
#'for(i in seq(0,100,by =20)){
#'  j <- i+20
#'
#'  list[[j/20]] <- SpatialPolygonsDataFrame(polygons(spec_poly)[i:j],spec_poly@data[i:j,])
#'}
#'
#'test <- combine_spdf(list[1:2])
#'@export
combine_spdf <- function(input,
         binomial = 'binomial'){

  for(i in seq(1,length(input),by =2)){
    if(i < length(input)){
    x <- input[[i]]
    y <- input[[i+1]]

    #add columns binomial if user specifies a different column
    if(binomial == 'binomial'){}else{
      y$binomial <- paste0("y$",binomial)
      x$binomial <- paste0("x$",binomial)
    }

    print('updating first spdf')
    #Merge x and y together to get consistent columns (technicalklky not needed)
    xy_merge <- merge(x@data, y@data, by = 'binomial', all.x = TRUE)
    xy_merge <- xy_merge[!duplicated(xy_merge$binomial),]

    print('updating second spdf')

    #Merge y and x together to get consistent columns
    yx_merge <- merge(y@data, x@data, by = 'binomial', all.x = TRUE)
    yx_merge <- yx_merge[!duplicated(yx_merge$binomial),]

    #used to update data
    x@data <- xy_merge
    y@data <- yx_merge

    #ensure cRS is the same (need to update)
    proj4string(x) <- proj4string(y)

    print('combining datasets')

    #combine
    spec_poly<- rbind(x,y)
    if(i == 1){spec_fin <- spec_poly}else{
      spec_fin <- rbind(spec_fin,spec_poly)
    }} else {
      x <- input[[i]]
      if(binomial == 'binomial'){}else{x$binomial <- paste0("x$",binomial)}
      x@data <- x@data$binomial
      spec_fin <- rbind(spec_fin,x)}

  }

print('aggregating')
#aggregate by binomial
#spec_poly_agg <- raster::aggregate(spec_fin, by = 'binomial',dissolve = TRUE)

print('buffering')
#add a buffer of length 0 (solves geometry errors often caused by merging datasets)
#spec_poly_b <- rgeos::gBuffer(spec_poly_fin, byid = TRUE, width = 0)

return(spec_fin)
}



