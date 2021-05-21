#Functions are stored here

###Cleaning RawGBIF Data####

clean_gbif <- function(df, filter = NULL, resolveTaxonomy = FALSE, resolveID = NULL){

  require(tidyverse)
  require(CoordinateCleaner)
  require(taxize)
  #filter by Genus

  print(nrow(df))
  df <- df[,c("gbifID","basisOfRecord","decimalLatitude","decimalLongitude",
              "coordinateUncertaintyInMeters","coordinatePrecision",
              "scientificName",
              "class","order","family","genus","species")]

  if(!is.null(filter)){

    df <- filter(df, df$genus == filter)

    print(nrow(df))
  }

  #clean blank species
  df <- df %>%
    filter(!is.na(df$species))
  df <- df%>%
    filter(df$species!= ",")
  df <- df%>%
    filter(df$species != "")
  df <- df %>%
    filter(df$species != " ")

  #remove blank coordinates
  df <- df %>%
    filter(.$decimalLatitude != ",")%>%
    filter(.$decimalLongitude != ",")

  #remove records such as fossils
  df <- filter(df, df$basisOfRecord == "HUMAN_OBSERVATION" |
                 df$basisOfRecord == 'OBSERVATION' |
                 df$basisOfRecord == 'PRESERVED_SPECIMEN')

  df <- filter(df, is.numeric(df$coordinateUncertaintyInMeters) <= 50000)

  names(df)[3:4] <- c('decimallatitude', 'decimallongitude')

  #clean using CoordinateCleaner
  df <- df %>%
    cc_val() %>%
    cc_equ() %>%
    cc_cap() %>%
    cc_cen() %>%
    cc_gbif() %>%
    cc_inst() %>%
    cc_sea() %>%
    cc_zero() %>%
    cc_dupl()

  return(df)
}

###CLeaning Data with KDE####
KDE_filter <- function(points, unispecies = TRUE, low_r = .25, up_r = .75){
  ###Clean using Kernel Density Estimates according to gomez 2018####
  require(spatstat)
  require(tidyverse)
  require(raster)

  if(unispecies == TRUE){
    #points_clean_kde = data.frame(x,y)

    win <- extent(matrix(c(points$decimallongitude,points$decimallatitude), nrow = nrow(points)))
    win <- data.frame(c(win@xmin, win@xmax),
                      c(win@ymin, win@ymax))

    win <- owin(c(win[1,1],win[2,1]),c(win[1,2],win[2,2]))

    #create ppm
    gbif_ppp <- ppp(points$decimallongitude,points$decimallatitude, window = win)

    #KDE
    gbif_kde <- density.ppp(gbif_ppp, at = "points")
    gbif_kde <- data.frame(gbif_kde, gbif_ppp$x,gbif_ppp$y)

    #Exclude Outliers
    q <- quantile(gbif_kde$gbif_kde, probs = c(low_r,up_r))
    iqr <- IQR(gbif_kde$gbif_kde)
    uq <- q[1] + iqr
    lq <- q[2] - iqr

    gbif_kde <- gbif_kde[gbif_kde$gbif_kde <= uq,]
    gbif_kde <- gbif_kde[gbif_kde$gbif_kde >= lq,]

  }else{

    points_clean_kde = data.frame("species" = as.character(),
                                  'x' = as.numeric(),
                                  'y' = as.numeric(),
                                  "KDE_status" = as.character())

    for (i in seq_along(levels(factor(points[,1])))){
      print(i)

      points_ <- subset(points, points[,1] == levels(factor(points[,1]))[i])

      if(nrow(points_)<=5){

        print(paste("Removing ", points_[1,1], " due to fewer than 5 records"))
        points_merge <- data.frame('species' = points_[1,1],
                                   'x' = points_$decimallongitude, 'y'= points_$decimallatitude,
                                   "KDE_status" = "Removed_tooFewRecords")

        points_clean_kde <- rbind(points_clean_kde,points_merge)
        next
      }else{

        win <- extent(matrix(c(points_$decimallongitude,points_$decimallatitude), nrow = nrow(points_)))
        win <- data.frame(c(win@xmin, win@xmax),
                          c(win@ymin, win@ymax))
        win <- owin(c(win[1,1],win[2,1]),c(win[1,2],win[2,2]))

        #create ppm
        gbif_ppp <- ppp(points_$decimallongitude,points_$decimallatitude, window = win)

        #KDE
        gbif_kde <- density.ppp(gbif_ppp, at = "points")
        gbif_kde <- data.frame(gbif_kde, gbif_ppp$x,gbif_ppp$y)

        #Exclude Outliers
        q <- quantile(gbif_kde$gbif_kde, probs = c(low_r,up_r))
        iqr <- IQR(gbif_kde$gbif_kde)
        uq <- q[2] + (iqr * 1.5)
        lq <- q[1] - (iqr * 1.5)

        gbif_kde <- gbif_kde[gbif_kde$gbif_kde <= uq,]
        gbif_kde <- gbif_kde[gbif_kde$gbif_kde >= lq,]

        nrow(gbif_kde)

        if(nrow(gbif_kde)<=1){

          print(paste(points_[1,1], " removed following filtering "))
          points_merge <- data.frame('species' = points_[1,1],
                                     'x' = points_$decimallongitude, 'y'= points_$decimallatitude,
                                     "KDE_status" = "Removed_KDE_RemovedAll")
          points_clean_kde <- rbind(points_clean_kde,points_merge)
          next
        }else{
          print(paste('Adding ',points_[1,1]))

          points_merge <- data.frame('species' = points_[1,1],
                                     'x' = gbif_kde$gbif_ppp.x, 'y'= gbif_kde$gbif_ppp.y,
                                     "KDE_status" = "Filter_Successful")
          points_clean_kde <- rbind(points_clean_kde,points_merge)

        }
      }
    }
  }
  return(points_clean_kde)
}

###Create circles around points and/or alpha polygons(NEED TO TEST)#####
# Points as a dataframe with columns species, x, and y
# map created from map function (used to crop the circles to a landmass)
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

###Create map####
mapcreate <- function(country_names){

  require(maps)
  require(maptools)

map <- maps::map('world',
                 regions = country_names, #c("Mexico",
                             #'Belize',
                             #'Guatemala',
                             #'Honduras',
                             #'El Salvador',
                             #'Costa Rica',
                             #'Panama',
                             #'Nicaragua'),
                 fill = TRUE,
                 interior = FALSE)


IDs <- sapply(strsplit(map$names,":"),function(x) x[1])

map_sp <- maptools::map2SpatialPolygons(map, IDs = IDs)%>%
  raster::aggregate(., dissolv = TRUE)

return(map_sp)

}

###Checks what species intersect?###
Iterate_intersect <- function(spdf_in, raster_size = numeric()){

  require('sf')
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

### Find optimal number of divisions for Allopatry ###
optimal_div <- function(spdf_in, min,max) {
  require('sf')
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


### Creates list of where species intersect with what chunks###
chunk_pres_ab <- function(spdf_in, raster_size = numeric(),species_spdf){

  require('sf')
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



  sp_inter <- sf::st_intersects(x_sf[],spdf_species[])
  return(list(inter,sp_inter))

}

