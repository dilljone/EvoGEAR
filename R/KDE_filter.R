#' Clean records based off kernel density estimate filter
#'
#'@param points a dataframe of records
#'
#'@export
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
