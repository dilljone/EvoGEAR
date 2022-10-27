#' Clean raw GBIF data
#'
#'@param df a dataframe of raw GBIF records to filter
#'@param filter An optional filter by binomial name
#'
#'@export
clean_gbif <- function(df, filter = NULL){

  require(CoordinateCleaner)
  #filter by Genus
  df <- GBIF_herps
  print(nrow(df))
  df <- df[,c("gbifID","decimalLatitude","decimalLongitude",
              "coordinateUncertaintyInMeters","coordinatePrecision",
              "scientificName", "basisOfRecord",
              "class","order","family","genus","species")]

  if(!is.null(filter)){

    df <- filter(df, df$genus == filter)

    print(nrow(df))
  }


  df <- df %>%
    filter(!is.na(species))%>%
    filter(species!= ",")%>%
    filter(species != "")%>%
    filter(species != " ")%>%
    filter(decimalLatitude != ",")%>%
    filter(decimalLongitude != ",")%>%
    filter(basisOfRecord == "HUMAN_OBSERVATION" |
             basisOfRecord == 'OBSERVATION' |
             basisOfRecord == 'PRESERVED_SPECIMEN')

  df <- filter(df, is.numeric(df$coordinateUncertaintyInMeters) <= 50000)

  names(df)[2:3] <- c('decimallatitude', 'decimallongitude')

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
