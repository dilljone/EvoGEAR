taxa <- unique(species_gbif_Iucn$x)

validateTaxa <- function(taxa){

#using taxize to clean the scientific names and merge

    syntaxa  <- taxize::synonyms(taxa, db = "eol", rows = 1)

    df_na <- ""

    df_tax <- data.frame(sub_tsn=character(),acc_name=character(),
                         acc_tsn=character(),acc_author=character(),
                         syn_author=character(),syn_name=character(),syn_tsn=character())

    df_syn <- data.frame(sub_tsn =character(),acc_tsn=character(),
                         syn_author=character(),syn_name=character(),syn_tsn=character())

    for(i in 1:length(syntaxa)){

      if(i==1){bool = 0}

      if(is.na(syntaxa[i])){

        #vector of species not found in ITIS
        df_na[length(df_na)+1] <- names(syntaxa[1])

      }else{

        temp <- data.frame(syntaxa[i])

        if(ncol(temp) == 7){
          print('here')
          if(bool == 0){

            #needed to make temp and change names because each df has unique column headers
            df_tax_temp = temp
            names(df_tax_temp) = c('sub_tsn','acc_name','acc_tsn','acc_author','syn_author','syn_name','syn_tsn')
          }

          else{df_tax <- rbind(df_tax,setNames(temp,names(df_tax_temp)))}

          bool = bool + 1

        }else{
          #logic for species found. Currently have logic for correct name but with synonyms. No code present for
          #correct name but no synonyms (returns 0col/0row dataframe)

          #df_syn_temp = temp
          #names(df_syn_temp) = c('sub_tsn','acc_tsn','syn_author','syn_name','syn_tsn')
          #df_syn <- rbind(df_syn,setNames(temp,names(df_syn_temp)))
        }
      }
    }
    names(df_tax)[6] = 'species'



  return(list(df_syn,df_tax,df_na))

}

?taxize::synonyms
test <- validateTaxa(phylo_species[1:1000,])


test<- taxize::synonyms(phylo_species[1:1000,], db = "eol")

save(out_phylo_0_1k,out_phylo_1_2k,out_phylo_2_3k,out_phylo_3_4k,out_phylo_4_5k,out_phylo_6_7k,out_phylo_7_8k,o)
