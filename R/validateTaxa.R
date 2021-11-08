validate_taxa <- function(input = as.list()){
  require(tidyverse)

  temp <- data.frame("binomial" = as.character(),
                     "synonym" = as.character(),
                     'source' = as.character())

  for(i in 1:length(input)){

    temp <- rbind(temp, data.frame(binomial = input[[i]]$binomial,
                                   synonym = input[[i]]$synonym,
                                   source = rep(names(input[i]),nrow(input[[i]]))))
  }

  temp <- unique(temp)

  output <- data.frame(matrix(nrow = length(unique(temp$binomial)), ncol = length(input)+1))
  names(output)[1] <- "Input_name"
  names(output)[2:ncol(output)] <- names(input)

  output$Input_name <- unique(temp$binomial)

  for(i in 1:nrow(output)){
    for(j in 1:length(input)){
      if(output[i,1] %in% input[[j]]$binomial){
        output[i,j+1] = 1
      }else(output[i,j+1] = 0)
    }
  }

  temp <- temp[which(is.na(temp$synonym)| temp$synonym != ""),]

  return(list(output,temp))
}

taxa_discrepency <- function(input = as.data.frame()){

  #requires input list where first column is accepted name, second column is synonym, third column is source
  #need to add a third column check skip

  require(tidyverse)

  temp = input
  rownames(temp) <- 1:nrow(temp)
  #create blank data frame
  matches <- data.frame(binomial= as.character(),
                        synonym = as.character(),
                        source = as.character())

  #iterate through input list
  for(i in 1:nrow(temp)){

    row_num <- NA
    row_num <- match(temp[i,2],temp$binomial)

    if(is.na(row_num)){}else{

      matches <- rbind(matches, temp[i,])
      matches <- rbind(matches, temp[row_num,])
    }
  }

  matches$note <- rep("",nrow(matches))
  matches <- matches %>%
    distinct(binomial,synonym,.keep_all = TRUE)

  for(i in 1:nrow(matches)){
    row_num <- NA
    row_num <- match(matches[i,2],matches$binomial)
    if(is.na(row_num)){  }else{
      if(matches[i,3] == matches[row_num,3]){
        matches[i,4] = "check"
      }else{
        matches[i,4] = "Synonym discrepency"
      }}

    if(matches[i,4] == ""){matches[i,4] = "Binomial discrpency"}

    if(nrow(matches[which(matches$binomial == matches[i,2]),])>1){
      matches[i,4] = "Remove. Cannot discern accepted name."
    }

  }

  return(matches)
}
