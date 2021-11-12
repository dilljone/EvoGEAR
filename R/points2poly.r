
#'Create circles around points and/or alpha polygons(NEED TO TEST)
#'
#'Points as a dataframe with columns species, x, and y
#'map created from map function (used to crop the circles to a landmass)
#'
#'

#gbif_raw <- read.delim("data/gbif.csv", header = TRUE, sep="\t")
#save(gbif_raw, file = "data/gbif.Rdata")

#load_all()


df <- load('data/gbif.Rdata')



  points2Poly <- function(points,x = 'x', y = 'y', species = 'species', map,
                        alpha = FALSE, a = 1,
                        type = 'All', buffer_m = 1000){
  require(raster)
  require(alphahull)
  require(dismo)
  require(rangeBuilder)

  spdf <- points

  #create spatialpointsdataframe

  coordinates(spdf) <- ~x+y
  projection(spdf) <- projection(map)
  spdf@data$species <- points$species

  #filtering records that are not on map (some species may be removed)

  spdf <- raster::intersect(spdf,map)

  #make circles around each point. D is in meters. Convert to polygon

  poly_list <- list()
  sp_list <- unique(spdf$species)

  for( i in 1:length(sp_list)){
    print(paste("Now doing", sp_list[i], i,"/",length(sp_list)))
    pol <- polygons(dismo::circles(spdf[which(spdf@data$species == sp_list[i]),],
                                   d=buffer_m,lonlat = TRUE))
    projection(pol) <- projection(map)
    #pol <- raster::intersect(pol,map)
    pol <- raster::union(pol)

    if(alpha == TRUE){
      poly_coords <- matrix(ncol = 2)

      for( i in 1:length((pol[1]@polygons[1][[1]]@Polygons))){
        poly_coords <- rbind(poly_coords, coordinates(pol[1]@polygons[1][[1]]@Polygons[i][[1]]))
      }

      poly_coords <- na.omit(unique(poly_coords))
      pol <- getDynamicAlphaHull(poly_coords,.90,1,buff = 1000, 3, clipToCoast = 'no')
      pol <- intersect(pol[[1]], map)
      poly_list[i] <- pol
      poly_list[[i]]@polygons[[1]]@ID = as.character(i)

    }
    else{}
    poly_list[i] <- pol
    poly_list[[i]]@polygons[[1]]@ID = as.character(i)
  }

  joined = SpatialPolygons(lapply(poly_list, function(x){x@polygons[[1]]}))
  jdata = SpatialPolygonsDataFrame(Sr=joined, data=data.frame(ID = 1:length(sp_list),species = sp_list),FALSE)

  return(jdata)
}
