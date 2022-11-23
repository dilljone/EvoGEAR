#'vis_acr_shiny is meant to launch a shiny application that helps you visualize your Ancestral State Reconstructions
#'
#'@param phylo is an object of class phylo with data related to states and their corresponding nodes
#'@param sf is an object of class sf that contains the states for visualization
#'@param state_list is a list of states and their frequencies. It should be named according to your nodes, and names must be just numbers (e.g. not node1, node2, node3 etc. Each node, needs an entry in the state list. Each dataframe in the list, should only have 2 columns, named state and freq in that order
#'@param most_likely is a TRUE/FALSE value. TRUE indicates that only the most likely freq is shown. FALSE indicates that the probability is instead used to visualize
#'
#'@details For the phylo object, each node should have a corresponding list that represents the states and their frequencies.
#' A common example output would be similar to that one provided by RASP. This list, should contain a single dataframe, with 2 columns.
#' The first column is the state, and the second column is the frequency.
#'
#'
#'
#'@export
#'

#where I left off, need to supply the list. Need to visualize with pie charts n shit. And do a better job so that you can visualize your stuff

set.seed(123456)
phylo <- ape::rtree(n = 20)
phylo$node =1:39

states = list()

for(i in 1:39){

  states[i] <- list(data.frame(state = sample(c("A",'B','C','D','AB','AC','AD','BC','BD','CD'),
                                                 5, replace = FALSE),
                                  freq = sample(c(1,2,3,4,5),5, replace = FALSE)/15))

  names(states)[i] <- paste0(i)

}

#Make the SF test
n1 = 2

data1 <- data.frame(x = rep(1:n1, each = n1),    # Create data frame for raster
                    y = rep(1:n1, n1),
                    value = runif(n1^2))

raster::rasterFromXYZ(data1)%>%
  as(.,"SpatialPolygonsDataFrame")%>%
  sf::st_as_sf() -> regions

regions$state <- c('A','B','C','D')

state_df <- ASR_list_2_df(states, ML = ML)

vis_acr_shiny(phylo, regions, state_df)

vis_acr_shiny <- function(phylo, sf, state_df){

  #Load required packages

  require(shiny)
  require(sf)
  require(tidyverse)
  require(plotly)
  require(ggtree)

  #blank states
  states_filt = ""

  #call ASR_list_2_df to convert a list of dataframes to a single dataframe

  gg_phylo<-ggtree(phylo)%<+% state_df

# UI ####
  ui <- fluidPage(

    # Application title
    title("ASE Visualizer"),
    hr(),
    selectInput('node',"Node Number:",gg_phylo$data$node),
    hr(),
    fluidRow(column(7,
                    plotlyOutput("phylo")),

             column(5,
                    plotOutput('state'))
    )
  )
# Server ####
  server <- function(input, output) {

    output$phylo <- renderPlotly({

      p2 <- gg_phylo +
        geom_point(aes(x = x,
                       y = y,
                       colour = ML_state,
                       label = node))

      plotly::ggplotly(p2, tooltip = c('node','ML_state'))

    })

    states <- reactive({

      states_filt <- gg_phylo$data%>%
        filter(node == as.numeric(input$node))%>% pull(ML_state)%>%
        str_split(.,"")%>%unlist()

      reg <- st_drop_geometry(sf)

      for(i in 1:nrow(reg)){

        if(reg[i,"state"] %in% states_filt){

          reg[i,"col"] = "Present"

        }else{
          reg[i,"col"] = "Not Present"}
      }

      sf$col = reg$col


      sf%>%
        #filter(col == "Present")%>%
        st_make_grid(.,cellsize = .5, square = FALSE) -> sf_pres

      #THIS IS ALL FUCKED UP. NEED TO FIX

      p <- ggplot(data = sf)+
        geom_sf(data = sf_pres, aes(fill = num))+
        geom_sf_label(aes(label = state))+
        scale_fill_gradientn(name = "regions", colours = terrain.colors(3))

      p

      return(p)
    })


    output$state <- renderPlot({
      states()

    })


  }

  # Run the application
  shinyApp(ui = ui, server = server)
}
