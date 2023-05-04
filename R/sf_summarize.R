#'This script will take in some SF object, have some grouping variable input, and output a series
#'of summary statistics that are relevant to biogeographical analyses
#'
#' @sf is the SF object in question.
#' @group is the name of the column that is used to group along the sf
#'
#'
#'
#'
#'


cluster_df <- data.frame(
  cluster = as.numeric(unique(hydro_phylo$cluster)),
  n_taxa = NA,
  n_edge = NA,
  mean_edge = NA,
  n_lc = NA,
  n_nt = NA,
  n_vu = NA,
  n_en = NA,
  n_cr = NA,
  n_hydro = NA,
  area_km2 = NA)

pull_clust <- function(cluster_df, sf, pres_ab){


#Build ot the database
for(i in 1:nrow(cluster_df)){

  c <- cluster_df[i,'cluster']

  #Added a tryCatch skip to account for regions with a low number of units (e.g. 1)

  skip <- FALSE

  x <-  tryCatch({pull_clust(c, hydro_phylo,pres_ab)}, error = function(e){skip <<- TRUE})

  if(skip){next}

  cluster_df[i,'n_taxa'] <- nrow(x)%>%as.numeric()

  cluster_df[i,'n_edge'] <- threat_levels%>%
    filter(edge_species == "EDGE species")%>%
    filter(taxon %in% rownames(x))%>%
    nrow() %>% as.numeric()

  cluster_df[i,'mean_edge']<- threat_levels%>%
    filter(taxon %in% rownames(x))%>%
    summarize(mean(EDGE_scores))%>%
    as.numeric()

  cluster_df[i,'n_lc'] <- threat_levels%>%
    filter(taxon %in% rownames(x))%>%
    filter(Threat_Code == "LC")%>%
    nrow() %>% as.numeric()

  cluster_df[i,'n_nt'] <- threat_levels%>%
    filter(taxon %in% rownames(x))%>%
    filter(Threat_Code == "NT")%>%
    nrow() %>% as.numeric()

  cluster_df[i,'n_vu'] <- threat_levels%>%
    filter(taxon %in% rownames(x))%>%
    filter(Threat_Code == "VU")%>%
    nrow() %>% as.numeric()

  cluster_df[i,'n_en'] <- threat_levels%>%
    filter(taxon %in% rownames(x))%>%
    filter(Threat_Code == "EN")%>%
    nrow() %>% as.numeric()

  cluster_df[i,'n_cr'] <- threat_levels%>%
    filter(taxon %in% rownames(x))%>%
    filter(Threat_Code == "CR")%>%
    nrow() %>% as.numeric()

  cluster_df[i,'n_hydro'] <- ncol(x)%>%as.numeric()

  cluster_df[i,'area_km2'] <- hydrobasin%>%
    filter(HYBAS_ID %in% colnames(x))%>%
    pull(SUB_AREA)%>%sum(.)

}


r_out <- data.frame(cluster = as.numeric(),
                    bio1_mean = as.numeric(),
                    bio1_median= as.numeric(),
                    bio1_sd = as.numeric(),
                    bio1_min = as.numeric(),
                    bio1_max = as.numeric(),
                    bio4_mean = as.numeric(),
                    bio4_median= as.numeric(),
                    bio4_sd = as.numeric(),
                    bio4_min = as.numeric(),
                    bio4_max = as.numeric(),
                    bio12_mean = as.numeric(),
                    bio12_median= as.numeric(),
                    bio12_sd = as.numeric(),
                    bio12_min = as.numeric(),
                    bio12_max = as.numeric())

for(i in 1:nrow(cluster_df)){

  c <- cluster_df[i,'cluster']
  #Added a tryCatch skip to account for regions with a low number of units (e.g. 1)

  skip <- FALSE
  print(c)
  x <-  tryCatch({pull_clust(c, hydro_phylo,pres_ab)}, error = function(e){skip <<- TRUE})

  if(skip){
    vect = c(c,NA,NA,NA,NA,NA)
    r_out <- rbind(r_out, vect)
    next}

  hydro_x <- hydro_simp%>%
    filter(HYBAS_ID %in% colnames(x))%>%
    mutate(ID = 1)%>%
    evogear::sf_collapse(.,"ID")%>%
    extract(bio,.)%>%
    do.call(rbind,.)%>%
    as.data.frame()%>%
    dplyr::select(bio1,bio4,bio12)

  foo <- function(x){
    out <- c(c, mean(x),median(x),sd(x),min(x),max(x))
    return(out)
  }

  vect <- vector()

  for(j in 1:(ncol(hydro_x))){

    vect <- c(vect, foo(hydro_x[,j]))

  }

  r_out <- rbind(r_out, vect)

}

names(r_out) = c('cluster','bio1_mean','bio1_median','bio1_sd','bio1_min','bio1_max',
                 'cluster_','bio4_mean','bio4_median','bio4_sd','bio4_min','bio4_max',
                 'cluster__','bio12_mean','bio12_median','bio12_sd','bio12_min','bio12_max')
