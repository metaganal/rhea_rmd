library(tidyverse)

run_date <- gsub(Sys.Date(), pattern = '-', replacement = '', fixed = T)


ar <- commandArgs(trailingOnly = TRUE)


rhea_generic_raw <- read_tsv(ar[[1]])


rhea_generic_chebi <- rhea_generic_raw %>%
  mutate(chebi = gsub('.*\\/', '', chebi)) %>%
  select(compoundAc, chebi) %>%
  distinct_all()

write_tsv(rhea_generic_chebi,
          path = paste0('data/', run_date, '_rhea_generic_to_chebi.tsv'))
