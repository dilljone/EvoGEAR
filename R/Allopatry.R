#redoing my allopatry script
#Steps:
# Goal: Create an adjacency matrix of degree of overlap
# eq: area of intersect/area of smaller Fitzpatrick and Tulleri 2006

#goal: Create a matrix of Euclidean distance between points
# for each occurence i of species x
#eq1: 0(xi) = w(xi)/b(xi) where w is ueclidian distance to nearest, and and b is distance to nearest heterospecif
#for each species
#eq2: p = n(0>1)/N(0)
#where n(0>1) is number of 0 values in which the distance to conspecific is greater than hetero and
#N(0) is total number of 0 values for species

#eq3: Thus, Overlap of each species is 0xy = (p(x) + p(y))/2
#0-1 bounded, where 0 is no overlap, and .5 is overlap is everything is random


#First Step. Adjacency matrix of degree overlap
# General Idea is to iterate over the SPDF and create degree overlap matrix for all species
#
# for (i in 1:10){
#   for(j in 1:10){
#
#     test <- raster::intersect(spec_poly[i,],spec_poly[j,])
#     plot(spec_poly[spdf$binomial == "Agalychnis_callidryas",], col = 'red')
#     plot(spec_poly[spec_poly$binomial == "Bolitoglossa_mexicana",], col = 'blue', add = TRUE)
#     plot(test, col = 'purple', add = TRUE)
#   }
# }
#
#
# test_ <- sf::st_intersects(spdf[spdf$binomial == "Agalychnis_callidryas",],
#                            spdf[spdf$binomial == "Bolitoglossa_mexicana",])
# test_
# test_ <- raster::intersect(spec_poly[spec_poly$binomial == "Agalychnis_callidryas",],
#                            spec_poly[spec_poly$binomial == "Bolitoglossa_mexicana",])
#
#
# rownames(inter_matrix[,1]) <- spec_poly[1,]$binomial
# intersect_mat <- as.matrix()
# inter_matrix <- data.frame()
# colnames(inter_matrix)[1]  <- "test"
# inter <- raster::crop(spec_poly[1,],spec_poly[1,])
#
# test <- as(spec_poly, "sf")
# (test[4,]$binomial)
# list(poly_test$binomial,poly_test$binomial)
# rownames(inter_matrix)[1]
#
# spdf <- as(poly_test,"sf")
# sf::st_crs(spdf) <- "+init=epsg:4326"
#
# area_over <- sum(sf::st_area(suppressWarnings(sf::st_intersection(spdf[1507,1], spdf[1517,1]))))
# plot(spdf[1517,3])
# lwgeom::lwgeom_make_valid(spdf[1507,])
# sf::st_
#
# x <- raster::raster(ncol = 5, nrow = 5, ext = raster::extent(poly_test), vals = 1)
# plot(x)
# plot(poly_test[5,], add = TRUE)
#
#
# inter <- sf::st_intersects(x_sf[,], spdf[])
#
#
# x_sf <- as(x, "SpatialPolygonsDataFrame")%>%
#   st_as_sf()


Iterate_intersect <- function(spdf_in, raster_size = numeric()){

  # require(sf)
  # require(tidyverse)
  output_matrix <- matrix(nrow = 1, ncol = 1, dimnames = list("row","Column"))


  #Create A raster that divides the data into easier chunks based on raster size
  #This is a square raster
  x <- raster::raster(ncol = raster_size,
                      nrow = raster_size,
                      ext = raster::extent(spdf_in),
                      vals = 1)
  #we then convert that into a sf object
  x_sf <- as(x, "SpatialPolygonsDataFrame")%>%
    st_as_sf()


  #convert SPDF to sf
  spdf_raw <- as(spdf_in,"sf")
  sf::st_crs(spdf_raw) <- "+init=epsg:4326"
  sf::st_crs(x_sf) <- "+init=epsg:4326"
  #then we get which polygons are present in each spdf. this is a binary function
  #that will output a list, where each item in the list is a vector containing the rows that
  # from your sf object that are present in each grid cell
  inter <- sf::st_intersects(x_sf[], spdf_raw[])

  #Next we iterate through the list creating subsets of our spdf object for each grid cell
  #then perform the matrix calculations. Note that this process sitll contains the geometry that
  #goes outside the cells, but only contains geometries that intersect into each cell
  length(inter)

  #now so we dont do the same comparisons twice we will get a taxon list
  taxon <- spdf_raw$binomial
  taxon$ID <- rownames(taxon)


  for(k in 1:length(inter)){
    print(paste("working on section ", k))
    spdf <- spdf_raw[inter[[k]],]

    #intermatrix is intermediary output object
    inter_matrix <- matrix(nrow = nrow(spdf), ncol = nrow(spdf),
                           dimnames = list(spdf$binomial,spdf$binomial))


    for(i in 1:(nrow(spdf)-1)){
      print(paste(i,"/",nrow(spdf)))
      for(j in (i+1):nrow(spdf)){
        print(paste(j,"/",nrow(spdf)))
        print(paste("now doing ",k,": ", i," x ",j))
        if(sf::st_is_valid(spdf[i,1])){
          if(sf::st_is_valid(spdf[j,1])){
            area_over <- sum(sf::st_area(suppressWarnings(sf::st_intersection(spdf[i,1], spdf[j,1]))))
            area_smaller <- min(sum(st_area(spdf[i,1])),sum(st_area(spdf[j,1])))

            #area_larger <- max(st_area(spdf[i,1]),st_area(spdf[j,1]))}

            if(length(area_over)==0){
              area_over <- 0
            }else{

            }
            print(area_over)
            inter_matrix[i,j] <- area_over/area_smaller
            rownames(inter_matrix)[i] <- spdf[i,]$binomial
            colnames(inter_matrix)[j] <- spdf[j,]$binomial
          }else{
            print(paste('failed on', spdf[j,]$binomial))
            inter_matrix[i,j] <- paste('failed on', spdf[j,]$binomial)
            rownames(inter_matrix)[i] <- spdf[i,]$binomial
            colnames(inter_matrix)[j] <- spdf[j,]$binomial
          }
        }else{
          print(paste('failed on', spdf[i,]$binomial))
          inter_matrix[i,j] <- paste('failed on', spdf[i,]$binomial)
          rownames(inter_matrix)[i] <- spdf[i,]$binomial
          colnames(inter_matrix)[j] <- spdf[j,]$binomial
        }
      }
    }

    output_matrix <- merge(output_matrix,inter_matrix, all = TRUE)
    rownames(output_matrix)<- colnames(output_matrix)
  }

  return(output_matrix)
}


optimal_div <- function(spdf_in, min,max) {
  # require('sf')
  # require(tidyverse)
  output_mat <- matrix(ncol = 2, nrow = 1)
  #convert SPDF to sf
  spdf_raw <- as(spdf_in,"sf")
  taxon <- spdf_in$binomial


  for (i in min:max){
    #Create A raster that divides the data into easier chunks based on raster size
    #This is a square raster

    check_mat <- matrix(ncol = length(taxon), nrow = length(taxon),
                        dimnames = list(taxon,taxon))

    x <- raster::raster(ncol = i,
                        nrow = i,
                        ext = raster::extent(spdf_in),
                        vals = 1)
    #we then convert that into a sf object
    x_sf <- as(x, "SpatialPolygonsDataFrame")%>%
      st_as_sf()

    sf::st_crs(spdf_raw) <- "+init=epsg:4326"
    sf::st_crs(x_sf) <- "+init=epsg:4326"
    #then we get which polygons are present in each spdf. this is a binary function
    #that will output a list, where each item in the list is a vector containing the rows that
    # from your sf object that are present in each grid cell
    inter <- sf::st_intersects(x_sf[], spdf_raw[])

    totalcalc <- 0

    for (j in 1:length(inter)){
      n_calc = 0
      if (length(inter[[j]]>0)){
        for (k in 1:length(inter[[j]])){
          for (l in k:length(inter[[j]])){
            if (is.na(check_mat[inter[[j]][k],inter[[j]][l]])){
              check_mat[inter[[j]][k],inter[[j]][l]] <- 1
              n_calc = n_calc + 1
            }
          }
        }
        totalcalc <- totalcalc + n_calc
      }else{
        totalcalc = totalcalc
      }
    }

    output_mat <- rbind(output_mat, c(i,totalcalc))

  }
  output_mat <- data.frame("Num_Div" = output_mat[-1,1],
                           "Num_calc" = output_mat[-1,2])
  return(list(output_mat,check_mat))
}

# opt_div <- optimal_div(agg,1,50)
#
#
# check <- opt_div[[2]]
# opt_div <- opt_div[[1]]
# plot(opt_div[,1], opt_div[,2],
#      main = "Optimal number of pairwise calculation divisions",
#      xlab = "Number of divisions",
#      ylab = "Number of calculations to perform")
#
# arrows(x1 = opt_div[opt_div[,2] == min(opt_div[,2]),1],
#        y1 = min(opt_div[,2]), col = 'blue',
#        x0 = opt_div[opt_div[,2] == min(opt_div[,2]),1]-1,
#        y0 = min(opt_div[,2])-1)
#
# ?arrows
# opt_div[opt_div[,2] == min(opt_div[,2]),1]
# min(opt_div[,2])
# out <- Iterate_intersect(agg,5)
#
# poly_test_sub <- subset(poly_test, poly_test@data$id_no %in% unique(poly_test@data[,c(1,2,7)]$binomial)$id_no)
#
# poly_test <- poly_test[order(poly_test$binomial,poly_test$year_),]
#
# agg <- raster::aggregate(poly_test,by = c('id_no','binomial','year_'), dissolve = TRUE)
# agg <- agg[order(agg$binomial,agg$year_, decreasing = TRUE),]
# agg@data$dupe <- duplicated(agg$binomial)
# agg <- agg[agg$dupe == FALSE,]
#
# merge()
#
#
#
#
#
#
#
# length(poly_test)
#
# for(j in 1:nrow()){
#   for(i in 1:nrow(spec_poly)){
#     print(c(i,":",j))
#     inter <- raster::crop(spec_poly[i,],spec_poly[j,])
#
#     if(is.null(inter)){
#       area <- 0
#     }else{
#       area <- sum(raster::area(inter))
#     }
#     print(area)
#     inter_matrix[i,j] <-  area
#     rownames(inter_matrix)[i] <- spec_poly[i,]$binomial
#     colnames(inter_matrix)[j] <- spec_poly[j,]$binomial
#     area ==0
#   }
# }
#
# rownames(inter_matrix)
#
# poly_test <- spec_poly_AS[spec_poly_AS$binomial %in% phylo_AS_raw$tip.label,]
#
# testing <- sf::st_intersection(test[1,1], test[1,1])
# sf::st_area(testing)
# apply(spec_poly[[]],1, function(x) print(x))
#
# test <- as(spec_poly, "sf")
# plot(test[2,1])
#
# test <- sf::st_as_sf(get(spec_poly[[1]]))
# colname
