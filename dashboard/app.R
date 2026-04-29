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
thematic::thematic_shiny(font = "auto")
# setup fonts

# Define UI for application that draws a barplot
ui <- fluidPage(
  
  # Application title
  titlePanel("Enteric Virus Reads"),
  
  #Everyone's favorite - styling! bslib themes
  theme = bs_theme(),
  

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
                  choices = c("scientific_name", "sample_srx", "year"))
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Bar Plot", plotOutput("barPlot")),
        tabPanel("Organism Frequency", plotOutput("org_freq")),
        tabPanel("Organism Percent", plotOutput("org_percent")),
        tabPanel("Organism Location", plotOutput("org_locale")),
        tabPanel("Map", leafletOutput("map"))
      )))
    
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  #themer
  bs_themer()


  # Load data - pull most recent date
  most_recent_kraken <- list.files("../results/summary/*", full.names = TRUE) %>% 
  sort(decreasing = TRUE) %>% 
  first()
  
  kraken_summary <- read.csv("../results/summary/kraken_summary-2026-04-28.csv")
  
  ww_metadata <- read.csv("../results/sra/sra_meta.tsv", sep = "\t")
  colnames(ww_metadata) <- c("sample_srx","library_strategy","library_selection","run_taxon_id","run_scientific_name", "collection_date", "location", "lat_lon")

  summary_data <- kraken_summary %>%
    left_join(ww_metadata, by = c("sample_srx" = "sample_srx")) %>%
    mutate(country = str_extract(location, "^[^:]+")) %>%
    mutate(place = str_extract(location, "(?<=:).*"))  
  
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
  
  
}

#save html
save_html(ui, file = "enteric_ww_dash.html")
# Run the application 
shinyApp(ui = ui, server = server)
