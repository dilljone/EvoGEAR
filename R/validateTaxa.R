validateTaxa <- function(taxa){


#using taxize to clean the scientific names and merge

    syntaxa  <- taxize::synonyms(taxa)

    for(i in 1:length(syntaxa)){

      if(i==1){bool = 0}

      if(is.na(syntaxa[i])){
        #Enter logic for NA here
      }else{
        temp <- data.frame(syntaxa[i])

        if(ncol(temp) == 7){
          print('here')
          if(bool == 0){
            df_tax = temp
            names(df_tax) = c('sub_tsn','acc_name','acc_tsn','acc_author','syn_author','syn_name','syn_tsn')
          }
          else{df_tax <- rbind(df_tax,setNames(temp,names(df_tax)))}
          bool = bool + 1
        }else{
          #enter logic here
        }
      }
    }
    names(df_tax)[6] = 'species'

    df <- plyr::join(df,df_tax[c(6,2,3),], by = 'species', type='left')

  return(df)

}
