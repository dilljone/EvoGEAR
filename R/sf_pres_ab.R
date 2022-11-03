#'sf_pres_ab will create a presence absence matrix given 2 sf objects
#'
#' @param x the first sf object
#' @param y the second sf object
#' @param x.unit the units of the x sf object. must be a column in the object
#' @param y.unit the units of the y sf object. must be a column in the object
#'
#'
#'
#'
#'@export

sf_pres_ab <- function(x,y,x.unit,y.unit){
  require(sf)
  require(tidyverse)
  #add in checks for if function is sf or not
  #add in checks that the units exist


  inter <-st_intersects(x,y)
  #outputs a list where each item in list is X
  #the row number of y is given within each item x

  mat <- matrix(nrow = nrow(x),ncol = nrow(y))
  rownames(mat) <- pull(x,x.unit) #name the rows based on user input
  colnames(mat) <- pull(y,y.unit) #name the cols based on user input

  for(i in 1:length(inter)){

    for(j in 1:length(inter[[i]])){

      mat[i,inter[[i]][j]] = 1
    }
  }
  mat[is.na(mat)] <- 0

  return(mat)
}
