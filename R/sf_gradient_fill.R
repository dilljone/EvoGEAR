#'sf_gradient_fill creates a gradient fill for use in the vis_acr_shiny application
#'
#'@param sf is the sf object
#'@param grid is the supplied grid over which the gradient will be created. Grids can be created with only a tiny bit of code
#'
#'@export

sf_gradient_fill <- function(sf, grid){


  st_combine(sf)%>%
    st_centroid()%>%
    st_coordinates%>%
    as_tibble() -> center


    df <- sf%>%
      st_centroid()%>%
      st_coordinates()%>%
      as_tibble()%>%
      mutate(z = sqrt((.$X - center$X)^2 + (.$Y - center$Y)^2))


for(i in 1:nrow(grid)){

  grid[i,]%>%
    st_centroid()%>%
    st_coordinates()%>%
    as_tibble() -> point

        distances <- sqrt((df$X - point$X)^2 + (df$Y - point$Y)^2)
        grid[i,"col"] <- sum(distances * df$z)
}

    return(grid)

}
