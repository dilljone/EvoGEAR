% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_spdf.R
\name{combine_spdf}
\alias{combine_spdf}
\title{Combine a list of spdfs together and group by binomialk name}
\usage{
combine_spdf(input, binomial = "binomial")
}
\arguments{
\item{input}{a list of lists containing each spdf}

\item{binomial}{Variable relating to the name of the column where binomials}
}
\description{
Combine a list of spdfs together and group by binomialk name
}
\examples{
   x <- pp
y <- IUCN
y <- y[y@data$legend == 'Extant (resident)',]
y <- raster::aggregate(y, by = 'binomial',dissolve = TRUE)
x$binomial <-  x$species

list <- list()
for(i in seq(0,100,by =20)){
 j <- i+20

 list[[j/20]] <- SpatialPolygonsDataFrame(polygons(spec_poly)[i:j],spec_poly@data[i:j,])
}

test <- combine_spdf(list[1:2])
}
