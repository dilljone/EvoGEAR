#' Clean raw GBIF data
#'
#'@param df a dataframe of raw GBIF records to filter
#'@param filter An optional filter by binomial name
#'
#'@export
clean_gbif <- function(df, filter = NULL){

require(CoordinateCleaner)
#filter by Genus

print(nrow(df))
df <- df[,c("gbifID","decimalLatitude","decimalLongitude",
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
#df <- filter(df, df$basisOfRecord == "HUMAN_OBSERVATION" |
 #              df$basisOfRecord == 'OBSERVATION' |
  #             df$basisOfRecord == 'PRESERVED_SPECIMEN')

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
