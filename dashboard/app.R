#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(htmltools)
library(tidyverse)

# Define UI for application that draws a barplot
ui <- fluidPage(
  
  # Application title
  titlePanel("Enteric Virus Reads"),
  
  #UI input - sliders (date), button, arranging. - V 
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      selectInput("y_var",
                  "Y:",
                  choices = c("count", "percent")),
      selectInput("x_var",
                  "Group By:",
                  choices = c("scientific_name", "sample_srx", "year")),
      selectInput("fill_var",
                  "Color By:",
                  choices = c("scientific_name", "sample_srx", "year")),
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("barPlot")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Load data - pull most recent date
  most_recent_kraken <- list.files("results/summary", full.names = TRUE) %>% sort(decreasing = TRUE) %>% first()
  kraken_summary <- read.csv(most_recent_kraken)
  
  output$barPlot <- renderPlot({
    
    # draw the histogram with the specified number of bins
    kraken_summary %>%
      # Drop unclassified from results
      mutate(year = year(as.Date(summary_date))) %>% #using summary date just to test, but this needs to be changed to the collection date from the metadata
      filter(scientific_name != "unclassified") %>%
      ggplot() +
      geom_bar(aes_string(x = input$x_var, y=input$y_var, fill = input$fill_var), stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs( )
  })
  
  #import data and transform - Script from Hanley
  # tabulated data - data table on a tab
  # ggplot - taxon breakdown by state AND National, and frequency
  # ggplot - time bound
  # need to add metadata - IMPORTANT - year is currently the year summary was made - needs to be changed to collection year from metadata
  
  
  
  
  # mapping - Leaflet map. V 
  # state shading? organism
  
  
  
}

#save html
save_html(ui, file = "enteric_ww_dash.html")
# Run the application 
shinyApp(ui = ui, server = server)
