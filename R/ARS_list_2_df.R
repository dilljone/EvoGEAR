#' ARS_list_2_df converts a list of Ancestral State Reconstructions (such as the output from RASP) to a dataframe
#' This is used internally in the to be used with the ASR visualizer function
#'
#'
#'@param list a list of dataframes where each dataframe contains the states and the frequency for each node. Each element in the list should be named with the node name (e.g. node1, node2, node3) state column which has the states as alphabetical characters and no spaces (e.g. A, BC, ADBG) , and a freq column that states the probability of each state
#'@param ML a TRUE/FALSE stateing if you want to include a column for most likely state for each node. FALSE by default
#'
#'
#'
#'
#'@export
#'
#'
#'


ASR_list_2_df <- function(list, ML = FALSE){
  require(tidyverse)

 list %>%
   do.call('rbind',.)%>%
   as.data.frame()%>%
   distinct(state)%>%
   arrange(state)%>%
             pull(state) -> state

 df <- data.frame(nodes = names(list))

 df[,2:(length(state)+1)] = 0

 names(df) <- c("node",state)

 for(i in 1:length(list)){

  for (j in 1:nrow(list[[i]])){

    df[i,list[[i]][j,'state']] <- list[[i]][j,'freq']

  }
if(ML == TRUE){

  df$ML_state <- ""

  for(i in 1:nrow(ML_state_df)){

    df%>%
      dplyr::slice(i)%>%
      select(-"node")%>%
      t()%>%
      as.data.frame()%>%
      arrange(desc(V1))%>%
      slice(1)%>%
      rownames()-> df[i,"ML_state"]}

  df%>%relocate(node, ML_state) -> df

}

 }
 return(df)
}
