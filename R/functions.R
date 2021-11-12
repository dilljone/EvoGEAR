#Functions are stored here

###Cleaning RawGBIF Data####

clean_gbif <- function(df, filter = NULL){

  require(tidyverse)
  require(CoordinateCleaner)
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
KDE_filter <- function(points, unispecies = TRUE, low_r = .25, up_r = .75, rec_min = 5){
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

      if(nrow(points_)<= rec_min){

        print(paste("Removing ", points_[1,1], " due to fewer than ", rec_min, " records"))
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

###Calculate the instability index from Mesquite
instability_index <- function(phy_list = as.list()){
  #Need to get a dataframe that is taxa by distance
  #the Index for each taxon in in the set of trees
  #that is summed across all distances between trees
  #The data matrix can be rows of each species and the columnes
  #can be the same species replicated for all trees
  require(ape)

  I_mat <- matrix()

  for(i in 1:length(phy_list)){
    j = i+1

    if (j>length(phy_list)){
      break()
    }
    #create Dix and Diy and reorder them for consistency
    #Need to figure out if rowname sort is required but lol fuck it for now
    Dix<- ape::cophenetic.phylo(phy_list[[i]])%>%
      as.data.frame(.)%>%
      select(order(colnames(.)))

    Diy <- ape:: cophenetic.phylo(phy_list[[j]])%>%
      as.data.frame(.)%>%
      select(order(colnames(.)))

    numer <- abs(Dix-Diy)
    denom <- (Dix-Diy)^2

    comp_mat <- (numer/denom) %>%
      rename_all(.,function(x) paste0(colnames(.),"_phylo_",j))

    I_mat <- cbind(I_mat,comp_mat)

  }

  index <- sapply(I_mat[,-1], function (x) as.numeric(gsub("NaN",0,x)))
  index <- rowSums(index)
  names(index) <- rownames(I_mat)

  return(index)

}

### Calculate the Consensus Fork Index
calc_CFI <- function(phy_list,
                     ranges = c(.5,.75,.95)){
  require(ape)

  out_list <- list()
  list_names <- vector()
  for(i in 1:length(ranges)){
    con_tree <- consensus(phy_list,p = ranges[i])
    boot <- prop.clades(con_tree,phy_list)
    CFI <- Nnode(con_tree)/(length(con_tree$tip.label)-2)
    out_list <- append(out_list,CFI)
    list_names <- append(list_names,paste0("CFI_",gsub("0.","",
                                                       as.character(ranges[i]))))
  }
  names(out_list) <- list_names
  return(out_list)
}


#Detect rogue taxa and highlight. Uses the custom function instability_index() and calc_CFI()

detect_rogue_taxa <- function(phy_list = as.list(), #Feed in a list of bootstrap trees
                              method = "percent", #need a series of options
                              percent_threshold = .10,
                              drop_threshold = .1){

  require(ape)
  require(tidyverse)

  CFI_full <-  calc_CFI(phy_list)
  insta_index <- instability_index(phy_list)

  switch(method,
         "percent" = {
           #percent removes tips based on percent_threshold
           #and their instability index scores
           #i.e. a percent_threshold of .1 would remove the 10% of worst scoring tips
           insta_index <-data.frame(index = sort(insta_index),
                                    tips = names(insta_index))
           num_end <- ceiling((1-percent_threshold) * nrow(insta_index))
           tips_kept <- insta_index[1:num_end,]
           tips_removed <- insta_index[num_end+1:nrow(insta_index),]
           tree_sub = list()

           for (i in 1:length(phy_list)){
             tree_sub[[i]] <- drop.tip(phy_list[[i]],tips_removed$tips)}
           return(list(insta_index,tips_kept, tips_removed,tree_sub,num_end))

         },

         "no_drop" = {
           #no drop drops no taxa and just returns the indeces
           return(insta_index)},

         "iterative" = {
           #iterative drops the worst performing taxa
           #then recalculates CFI
           #then selectively drops taxa based on drop_threshold
           #need to add a way to not go through every single taxa

           insta_index <-data.frame(index = sort(insta_index),
                                    tips = names(insta_index))

           tips_remove <- data.frame(taxa = as.character(),
                                     CFI_50 = as.numeric(),
                                     CFI_75 = as.numeric(),
                                     CFI_95 = as.numeric())
           for(i in 1:nrow(insta_index))
           {
             #sub_tree <- drop.tip(phy_list[[i]],insta_index[i,2])
             tree_sub <- list()

             for (j in 1:length(phy_list)){
               tree_sub[[j]] <- drop.tip(phy_list[[j]],insta_index[i,2])
             }

             CFI_sub <- calc_CFI(tree_sub)

             CFI_diff_50 <- CFI_sub[[1]] - CFI_full[[1]]
             CFI_diff_75 <- CFI_sub[[2]] - CFI_full[[2]]
             CFI_diff_95 <- CFI_sub[[3]] - CFI_full[[3]]

             if(CFI_diff_95 >= drop_threshold |
                CFI_diff_75 >= drop_threshold |
                CFI_diff_50 >= drop_threshold){

               tips_remove <- rbind(tips_remove,
                                    data.frame(insta_index[i,2],
                                               CFI_diff_50,
                                               CFI_diff_75,
                                               CFI_diff_95))

             }
           }
           return(list(tips_remove,CFI_full))
         }
  )

  return(list(insta_index,tips_kept, tips_removed,tree_sub,num_end))
}

iterate_genbank <- function(search_terms, by = 1000,
                            seq_length_max = 30000,
                            extra_term = ""){

  require(rentrez)
  require(tidyverse)

  #Custom Sequence Function that ensures last sequences
  # If there are 369 sequences and you go by 100s, then the last 69 would not be included
  seq0 <- function(from = 1, to = 1, by = 1, incLast = TRUE){
    out = do.call(what = seq, args = list(from, to, by))
    if (incLast & to%%by != 0){
      out = c(out, tail(out, 1) + by)
    }
    return(out)
  }

  #iterate through each search term in list
  for(i in 1:length(search_terms)){

    term <- paste0("('",search_terms[i],"'","[Organism] OR ", search_terms[i],
                   "[All Fields])"," AND (1 [SLEN] : ",seq_length_max,"[SLEN])",extra_term)

    search_res <- entrez_search(db="nuccore", term=term,use_history = TRUE)
    #printing to make it look nicer. Adding progress bar as well

    print(paste0("Now Quering: ",search_terms[i], ". ", search_res$count, " records found"))
    pb <- txtProgressBar(min = 0, max = search_res$count, style = 3, width = 30, char = "\U1F438")

    #Ifelse to fetch if few sequences
    if(search_res$count<by){

      fetch <- entrez_fetch(db = "nuccore", #id = search_res$ids[1:length(search_res$ids)],
                            rettype = "fasta", web_history = search_res$web_history)

      write(fetch, paste0(search_terms[i],"_FASTA.fasta"), sep = "\n")
      setTxtProgressBar(pb,search_res$count)

    }else{
      #else is for larger amount of sequences
      start = by
      fetch <- entrez_fetch(db = "nuccore", #id = search_res$ids[1:length(search_res$ids)],
                            rettype = "fasta", retmax = by,
                            web_history = search_res$web_history)

      write(fetch, paste0(search_terms[i],"_FASTA.fasta"), sep = "\n")

      setTxtProgressBar(pb,by)

      #starting to pull sequences every by

      for(j in seq0(from = start,to = search_res$count, by = by)){

        setTxtProgressBar(pb,j)

        if(j==by){next}else{

          fetch <- entrez_fetch(db = "nuccore", #id = search_res$ids[start:j],
                                rettype = "fasta", web_history = search_res$web_history,
                                retmax = by, retstart = start)

          write(fetch, paste0(search_terms[i],"_FASTA.fasta"),sep = "\n", append = TRUE)
          start = j
        }
      }
    }
    #finish off the ProgressBar
    setTxtProgressBar(pb,search_res$count)
    close(pb)
  }
}

sp_list_amnh <- function(data){
  library(tidyverse)
  library(stringr)
  library(magrittr)

  data_out <-
    data %>%
    mutate("ID" = 1:nrow(data)) %>%
    mutate('Raw' = trimws(data[,1]))%>%
    #select('trimmed')%>%
    mutate('level' = substr(.$Raw,1,str_locate(.$Raw,":")-1))%>%
    mutate('value' = (trimws(str_split(.$Raw,":",n = 2, TRUE)[,2])))%>%
    mutate('Order' = ifelse(level == 'Order',value,'-'),
           'Superfamily' = ifelse(level == 'Superfamily',value,'-'),
           'Family' = ifelse(level == 'Family',value,'-'),
           'Subfamily' = ifelse(level == 'Subfamily',value,'-'),
           'Genus' = ifelse(level == 'Genus',value,'-'),
           'Species' = ifelse(level == 'Species',value,'-'))%>%
    mutate('binomial' = ifelse(level == 'Species',
                               (substr(Species,1,
                                       (str_locate(Species,"^([^ ]+[ ]+[^ ]+)[ ]")[,2]))),Species))%>%
    mutate(binomial = str_trim(binomial))%>%
    select(ID,Raw, level, Order, Superfamily,Family,Subfamily,Genus,Species,binomial)
  return(list(data_out,unique(data_out$binomial)))
}
