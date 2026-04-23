
# Load packages (install if necessary)
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}


library("tidyverse")


# Read in Kraken report & add SRX id (should this be SRA?) as a column
kraken_reports_files <- list.files("test/kraken2", full.names = TRUE)

kraken_reports <- bind_rows(lapply(kraken_reports_files, function(x) {
  read.delim2(x, header = FALSE) %>%
    mutate(
      source_file = gsub("\\.kraken2.*", "", basename(x)),
      date_summary = Sys.Date()
      )
}))


colnames(kraken_reports) <- c("percent", "count", "count_taxon_exclusive", "rank_code", "ncbi_taxon_id", "scientific_name", "sample_srx", "summary_date")
  

# identify taxa of interest based on fasta files in reference folder
ref_tax_ids <- list.files("./ref_fastas") %>%
  str_extract("(?<=NC_)[^.|^_]+")


# Filter kraken_reports
kraken_comb <- filter(kraken_reports) %>%
  filter(ncbi_taxon_id %in% ref_tax_ids | rank_code == "U") %>% # Extract only enteric viruses of interest (and also unclassified count)
  mutate(scientific_name = gsub("^\\s+", "", scientific_name)) 

# Convert to wide format (1 row per samples)
kraken_comb_wide <- kraken_comb  %>%
  select(sample_srx, summary_date, scientific_name, percent) %>%
  pivot_wider(
    names_from = "scientific_name",
    values_from = "percent")

# save files
if (!dir.exists("./test/summary")) {
  dir.create("./test/summary", recursive = TRUE)
}

write.csv(kraken_comb, paste0("./test/summary/kraken_summary-", Sys.Date(), ".csv"), row.names = FALSE)

write.csv(kraken_comb_wide, paste0("./test/summary/kraken_summary_wide_format-", Sys.Date(), ".csv"), row.names = FALSE)
