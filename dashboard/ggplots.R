setwd("C:/Users/emily/Documents/LACPHL_work/aphl_hackathon")

library(ggplot2)
library(readr)

summaryfile <- read_csv("kraken_summary-2026-04-28.csv") ##initial file read in - no metadata attached

summary_oi <- summaryfile[summaryfile$scientific_name != "unclassified", ] ##subset to exclude unclassified

ggplot(summary_oi, aes(y=scientific_name, fill=scientific_name )) + 
  geom_bar( ) +
  scale_fill_hue(c = 40) +
  theme(legend.position="none") + labs(title = "Count of organisms found in wastewater", x="Number of organisms", y="Organism")  ##bar chart to number organisms by scientific name

ggplot(summary_oi, aes(x=percent, y=scientific_name)) +
  geom_line( color= "salmon", linetype=1, size = 5) +
  labs(title="Percentage of Organisms in Wastewater", x="Percentage") ##bar chart for percent of organisms by scientific name (without unclassified)

ggplot(summary_oi, aes(x = location)) +
  geom_bar(fill = "salmon") +
  labs(title = "Organisms by Location", y = "Number of Organisms")  ##bar chart by location - should work if column name in Toms sheet is location lol
