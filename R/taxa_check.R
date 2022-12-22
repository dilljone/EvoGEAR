#' Checks taxa given a named list of taxa
#' Use case would be to supply taxa from your phylogeny, occurence records, and known databases
#' Eventually will attempt to add in already known taxa databases
#' #need to peroperly add at a later date
#'
#'
#'
#'
#'
#'
#'
#'
#'


taxa_check <- function(list = character()){

  taxa <- unlist(list)%>%unique()
  out <- data.frame(taxa = taxa)

  for(i in 1:length(list)){

    out[,i+1] <- out$taxa %in% list[[i]]

  }

  colnames(out) <- c("taxa", names(list))
  return(out)
}
