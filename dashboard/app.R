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
library(ggplot2)
library(readr)
library(bslib)
library(leaflet)
library(DT)


thematic::thematic_shiny(font = "auto")
# setup fonts

# Define UI for application that draws a barplot
ui <- fluidPage(
  
  # Application title
  titlePanel("Enteric Virus Reads"),
  
  #Everyone's favorite - styling! bslib themes
  theme = bs_theme(),
  

  #UI input - sliders (date), button, arranging. - V 
  
  tabsetPanel(
    
    # barplot tab
    tabPanel("Bar Plot",
    # Set up sidebar
    sidebarLayout(
      sidebarPanel(
        selectInput("y_var",
                    "Summarize Average Read:",
                    choices = c("count", "percent")),
        selectInput("x_var",
                    "Group By:",
                    choices = c("scientific_name", "sample_srx")),
        selectInput("fill_var",
                    "Color By:",
                    choices = c("scientific_name", "sample_srx", "location", "year", "library_selection"))
                    ),
        mainPanel(plotOutput("barPlot"))
    )),
    
    # Timeplot
    tabPanel("Time Plot",
    sidebarLayout(
      sidebarPanel(
        selectInput("y_var2",
                    "Summarize Average Read:",
                    choices = c("count", "percent"))
      ),
      mainPanel(plotOutput("timePlot"))
    )),

    # Additional tabs
    tabPanel("Organism Frequency", plotOutput("org_freq")),
    tabPanel("Organism Percent", plotOutput("org_percent")),
    tabPanel("Organism Location", plotOutput("org_locale")),
    tabPanel("Map", leafletOutput("map")),
    tabPanel("Data", dataTableOutput("data"))

    )
)

# Define server logic
server <- function(input, output) {
  #themer
  bs_themer()


  # Load data

  ### Generate combined Kraken summary file
  try(source("../dashboard/summarize_kraken_output.R", chdir = TRUE))
  
  ### Load Kraken summary file based on most recent date in file name
  most_recent_kraken <- list.files("../results/summary/", full.names = TRUE) %>% 
  sort(decreasing = TRUE) %>% 
  first()
  
  kraken_summary <- read.csv("../results/summary/kraken_summary-2026-04-28.csv")
  
  ### Metadata from SRA data pull
  ww_metadata <- read.csv("../results/sra/sra_meta.tsv", sep = "\t")
  colnames(ww_metadata) <- c("sample_srx","library_strategy","library_selection","run_taxon_id","run_scientific_name", "collection_date", "location", "lat_lon")

  ### Combine Kraken & metadata
  summary_data <- kraken_summary %>%
    left_join(ww_metadata, by = "sample_srx") %>%
    mutate(
      country = str_extract(location, "^[^:]+"),
      place = str_extract(location, "(?<=:).*"),
      collection_date = as.Date(collection_date, format = "%Y-%m-%d"),
      year = as.character(year(as.Date(collection_date)))
    )
  

  # Generate barplot
  output$barPlot <- renderPlot({
    
    summary_data %>%
      # Drop unclassified from results
      filter(scientific_name != "unclassified") %>%
      ggplot() +
      stat_summary(fun = mean, geom = "col", position = "dodge",
      aes_string(x = input$x_var, y = input$y_var, fill = input$fill_var)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs( )
  })
  


# Time trends
  output$timePlot <- renderPlot({

    summary_data %>%
      # Drop unclassified from results
      filter(scientific_name != "unclassified") %>%
      ggplot(aes_string(x = "collection_date", y = input$y_var2, color = "location")) + 
      geom_point() +
      geom_line() +
      facet_wrap(~scientific_name, ncol = 2) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs( )
  })


  #import data and transform - Script from Hanley
  # tabulated data - data table on a tab
  # ggplot - taxon breakdown by state AND National, and frequency
  # ggplot - time bound

summaryfile <- summary_data 
##initial file read in - no metadata attached

summary_oi <- summaryfile[summaryfile$scientific_name != "unclassified", ] ##subset to exclude unclassified

output$org_freq <- renderPlot({
  ggplot(summary_oi, aes(y=scientific_name, fill=scientific_name )) + 
  geom_bar( ) +
  scale_fill_hue(c = 40) +
  theme(legend.position="none") + labs(title = "Count of organisms found in wastewater", x="Number of organisms", y="Organism")  
  ##bar chart to number organisms by scientific name
})
output$org_percent <- renderPlot({
  ggplot(summary_oi, aes(x=percent, y=scientific_name)) +
  geom_line( color= "salmon", linetype=1, size = 5) +
  labs(title="Percentage of Organisms in Wastewater", x="Percentage") 
  ##bar chart for percent of organisms by scientific name (without unclassified)
})
output$org_locale <- renderPlot({
  ggplot(summary_oi, aes(x = location)) +
  geom_bar(fill = "salmon") +
  labs(title = "Organisms by Location", y = "Number of Organisms") 
   ##bar chart by location - should work if column name in Toms sheet is location lol
})
  
  
  
  # mapping - Leaflet map. V 
  # state shading? organism
#output$map <- renderLeaflet({
#leaflet() %>%
#addTiles() %>% # Add default OpenStreetMap tiles
#setView(lng = 0.3458, lat = 50.6253, zoom = 3) # Set initial view
#})
summary_data_map <- summary_data %>%
  filter(scientific_name != "unclassified") %>%
  mutate(lat_lon = str_replace(lat_lon, "N", ", -")) %>%
  separate_wider_delim(lat_lon, names = c("lat", "lng"), delim = ",", too_few = "align_start") %>%
  mutate(lng = str_replace(lng, "W", "")) %>%
  mutate(lat = str_replace(lat, "- ", "-")) %>%
  mutate(lng = str_replace(lng, "- ", "-"))%>%
  mutate(lat = jitter(as.numeric(lat), factor = 0.2), lng = jitter(as.numeric(lng), factor = 0.2))

enteric_count <- summary_data_map %>%
  group_by(scientific_name, location) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

enteric_map <- summary_data_map %>%
  left_join(enteric_count, by = c("scientific_name", "location"))

output$map <- renderLeaflet({
    enteric_map %>%
      leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron,
        options = providerTileOptions(noWrap = TRUE)
      ) %>%
      addCircleMarkers(popup = ~scientific_name,
      color = ~colorFactor("viridis", summary_data_map$scientific_name)(scientific_name),
      fillOpacity = 0.4, 
      radius = ~count.y * 15, group = ~scientific_name, clusterOptions = markerClusterOptions())
       })
output$data <- 
  renderDataTable({datatable(summary_data)})   
  
}

#save html
save_html(ui, file = "enteric_ww_dash.html")
# Run the application 
shinyApp(ui = ui, server = server)
