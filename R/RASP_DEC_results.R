#' RASP DEC results
#'@param file is the file path to the RASP output dispersal file
#'
#'
#'
#'
#'@export
#'
#'

read_RASP_dispersal <- function(file){

  require(tidyverse)
  text <- read_table(skip_empty_rows = FALSE,
                     col_names = FALSE,
                     file)

  text$X1 %>% as.vector() -> text

  which(is.na(text)) -> na_vect

  c(0,na_vect) -> na_vect

  df <- data.frame(text = "")
  for(i in 1:length(na_vect)){

    j <- na_vect[i]

    val <- paste(text[(j+1):(j+9)], collapse= " ")
    df[i,] <- val

  }

  df <- df[-nrow(df),]%>%as.data.frame()

  df%>%separate(col = 1, sep = " ",
                into= c('Node','Event','Dispersal','Vicariance',
                        'Extinction','Event_','Path','Prob','Probability'))%>%
    select(Node,Dispersal,Vicariance,Extinction,Path,Probability)%>%
    mutate(Node = gsub("NODE","",Node),
           Node = gsub(":","",Node),
           Dispersal = gsub("Dispersal:","",Dispersal),
           Vicariance = gsub("Vicariance:","",Vicariance),
           Extinction = gsub("Extinction:","",Extinction),
           ML_state = str_extract(Path,"^[^-]*")) -> df

  return(df)

}


