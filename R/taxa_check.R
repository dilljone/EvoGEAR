#'This function will check if taxa are present according to a user supplied list
#'Can be used to check where taxa are being excluded, where there are no overlaps,
#'or to validate taxa across a list
#'
#'@param list a list of taxa from various datasets. Can be either dataframes or a list of vectors
#'@param var_name NULL by default. If supplied datasets are in dataframe format, you can specify the variable names
#'@param gbif_backbone Boolean. FALSE by default. Will query the GBIF backbone and supply a secondary list detailing
#'
#'This function returns a dataframe where each item in the list constitutes a single column
#'The rownames contain a unieue vector of taxa found in all datasets
#'
#'
#'@export




taxa_check <- function(list, var_name = NULL, gbif_backbone = F) {
  #List is a list of named datasets
  #var names is a vector of the variable name used in order of each  element of the list

  taxa_all <- vector()
  if(is.null(var_name)){
    for( i in 1:length(list)){
      taxa_all <- c(taxa_all, list[[i]])
    }
  }else{
    for( i in 1:length(list)){
      taxa_all <- c(taxa_all, list[[i]][,var_name])
    }
  }
  #Keep only unique
  taxa_all <- unique(taxa_all)%>%sort()

  out <- matrix(nrow = length(taxa_all), ncol = length(list))%>%as.data.frame()
  rownames(out) <- taxa_all
  colnames(out) <- names(list)
  for(i in 1:nrow(out)){
    for(j in 1:ncol(out)){
      if(taxa_all[i] %in% list[[j]]){
        out[i,j] = TRUE
      }else{ out[i,j] = FALSE}
    }
  }

  if(gbif_backbone == TRUE){

    require(rgbif)
    gbif <- rgbif::name_backbone_checklist(taxa_all)

    out$gbif_status <- gbif$status
    out$gbif_confidence <- gbif$confidence
    out$gbif_rank <- gbif$rank

    return(list(out,gbif))

  }else{
    return(out)

  }


}

taxa_all <- c("Rana vaillanti", "Lithobates sylvaticus", 'iguana iguana')
