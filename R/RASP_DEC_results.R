#' RASP DEC results
#'
#'
#'
#'
#'
#'
#'

text <- read_delim(delim ="\r", "F:/My Drive/OneDrive/Research/Phyloregion_SDSU/Phyloregionalization_MA/Outputs/RASP_Squamate_DEC_results.txt")[,-2]

t_pos <- which(text== "[TAXON]")
p_pos <- which(text== "[TREE]")
r_pos <- which(text== "[RESULT]")


taxa <- separate(data = text[t_pos+3:p_pos-2,1],
           col = colnames(text[t_pos+3:p_pos-2,1]),
           sep = '\\t',
           into = c("tip",'species','state'))%>%as.data.frame()

tree <- text[p_pos+1,]%>%pull()%>%ape::read.tree(text = .)

results <- text[r_pos+2:nrow(text),1]
