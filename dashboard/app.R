#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Before running this script for the frist time, it may be necissary to run in the terminal:
# conda install -c conda-forge r-base
# conda install -c conda-forge r-tidyverse r-shiny r-leaflet r-thematic r-DT

library(shiny)
library(htmltools)
library(tidyverse)
library(ggplot2)
library(readr)
library(bslib)
library(leaflet)
library(DT)
library(scales)


thematic::thematic_shiny(font = "auto")
# setup fonts




  # Load data

  ### Generate combined Kraken summary file
  try(source("../dashboard/summarize_kraken_output.R", chdir = TRUE))
  
  ### Load Kraken summary file based on most recent date in file name
  most_recent_kraken <- list.files("../results/summary/", full.names = TRUE, pattern = "kraken_summary-") %>% 
  sort(decreasing = TRUE) %>% 
  first()
  
  kraken_summary <- read.csv(most_recent_kraken)
  
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
  



# Define UI for application that draws a barplot
ui <- fluidPage(
  
  # Application title
  titlePanel("Enteric Virus Reads"),
  
  #Everyone's favorite - styling! bslib themes
  theme = bs_theme(),
  

  #UI input - sliders (date), button, arranging. - V 
 
  tabsetPanel(
    
    # primary tab
    tabPanel("Primary Summary",
      fluidRow(
        h3("Enteric Viruses in Metagenomic Wastewater Samples", align = "center"),
        column(width = 6, plotOutput("org_freq")),
        column(width = 6, plotOutput("org_percent")),
        column(width = 6, plotOutput("org_locale")),
        column(width = 6, plotOutput("unique_viruses"))
    )),

    # barplot tab
    tabPanel("Barplot",
    # Set up sidebar
    sidebarLayout(
      sidebarPanel(
        selectInput("y_var",
                    "Summarize Average Read:",
                    choices = c("count", "percentage_in_sample")),
        selectInput("x_var",
                    "Group By:",
                    choices = c("scientific_name", "sample_srx")),
        selectInput("fill_var",
                    "Color By:",
                    choices = c("scientific_name", "sample_srx", "location", "year", "library_selection"))
                    ),
        mainPanel(plotOutput("barPlot"))
    )),
    
    # Additional tabs
    tabPanel("Map", leafletOutput("map")),
    tabPanel("Data", dataTableOutput("data"))

    )
)

# Define server logic
server <- function(input, output) {
  #themer
  bs_themer()


  # Generate barplot
  output$barPlot <- renderPlot({
    
    summary_data %>%
      # Rename variables
      mutate(
        `percentage_in_sample` = percent
      ) %>%
      # Drop unclassified from results
      filter(scientific_name != "unclassified") %>%
      ggplot() +
      stat_summary(fun = mean, geom = "col", position = "dodge",
      aes_string(x = input$x_var, y = input$y_var, fill = input$fill_var)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs( )
  })
  

summary_oi <- summary_data[summary_data$scientific_name != "unclassified", ] ##subset to exclude unclassified

output$org_freq <- renderPlot({
  ggplot(summary_oi, aes(x=count, y=scientific_name, fill=scientific_name)) + 
  geom_boxplot( ) +
  geom_boxplot(color = "salmon") +
  scale_fill_hue(c = 40) + 
  theme(legend.position="none") +
  labs(title = "Count of organisms found in wastewater", x="Number of organisms (reads)", y="Organism")  
  ##bar chart to number organisms by scientific name
})
output$org_percent <- renderPlot({
  ggplot(summary_oi, aes(x=percent, y=scientific_name, fill=scientific_name)) +
  geom_boxplot(color = "salmon") +
  scale_fill_hue(c = 40) +
  scale_x_continuous(labels = label_number()) +
  theme(legend.position="none") +
  labs(title="Percentage of Organisms in Wastewater (axis is not to 100!)", x="Percentage in sample") 
  ##bar chart for percent of organisms by scientific name (without unclassified)
})
output$org_locale <- renderPlot({
  ggplot(summary_oi %>% mutate(location = ifelse(location == "", "unknown", location)), aes(x = location)) +
  geom_bar(fill = "salmon") +
  labs(title = "Sample Locations", y = "Number of Samples") 
})
output$unique_viruses <- renderPlot({
  summary_data %>%
  group_by(sample_srx) %>%
  reframe(unique_pathogens = n()-1) %>%
  ggplot(aes(x = unique_pathogens)) +
  geom_histogram(fill = "salmon", binwidth = 1) +
  scale_x_continuous(breaks = seq(0, 15, by = 1)) +
  labs(
    title = "Total Unique Enteric Viruses Identified per Sample",
    x = "unique enteric viruses in sample",
    y = "Total Samples") 
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
