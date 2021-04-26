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

  #create taxon list
  taxon <- as.character(df$species)%>%
    unique()
  
  if(resolveTaxonomy == TRUE){
    print("Now Resolving Taxonomy")
    #using taxize to clean the scientific names and merge
    new_tax <- taxize::gnr_resolve(taxon, best_match_only = TRUE, 
                                   canonical = TRUE, http = 'post', data_source_ids = as.numeric(resolveID))
    
    df <-merge(df,new_tax, by.x = 'species', by.y = 'user_supplied_name')
  }

  
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


sp_list <- as.list()
spdf <- spec_poly_done
cunt <- map_sp

sp_list <- as.character(factor(spdf$binomial))

crop_list <- list()

for (i in 1:length(sp_list)){

  sub <- subset(spdf,spdf$binomial == sp_list[i])
  rast <- raster::raster(ext = extent(sub), res = 2)
  cunt_crop <- raster::crop(cunt,rast)
  crop_list[[i]] <- cunt_crop
}
  
names(crop_list) <- sp_list
library(raster)

map <- maps::map(database="world", fill = TRUE, interior = TRUE)
IDs <- sapply(strsplit(map$names,":"),function(x) x[1])

map_sp <- maptools::map2SpatialPolygons(map, IDs = IDs)%>%
  raster::aggregate(., dissolv = TRUE)

library(tidyverse)
setwd("C:/Users/Admin/OneDrive - San Diego State University (SDSU.EDU)/Research/Phyloregion_SDSU/Phyloregionalization_MA/Outputs")
load("prephylo.data")

plot(crop_list$`Agkistrodon piscivorus`)
plot(subset(spec_poly_done,spec_poly_done$binomial == 'Agkistrodon piscivorus'),
     col = 'red', add = TRUE)
projection(spdf)
projection(map_sp)

extent(sub)

plot(map_sp)
extent(map_sp)

