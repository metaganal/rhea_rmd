# Set working directory to the git repository
setwd("~/rhea_rmd")


# Load package
library(tidyverse)

run_date <- gsub(Sys.Date(), pattern = "-", replacement = "", fixed = T)

# Parse command line arguments
ar <- commandArgs(trailingOnly = TRUE)


# Read RHEA table
dframe <- readr::read_tsv(ar[[1]])


# Extract IDs from table strings
rhea_table <- dframe %>%
  dplyr::mutate(reactionSide = gsub(".*[0-9]{5}_", "", reactionSide)) %>%
  dplyr::select(rheaid,
                reactionEquation,
                reactionSide,
                compoundID,
                compoundName) %>%
  dplyr::arrange(rheaid,
                 reactionSide)


# Summarise RHEA IDs variable into one row
rhea_reactions <- rhea_table %>%
  dplyr::group_by(rheaid, reactionSide) %>%
  dplyr::summarise(compound = paste(list(compoundID), collapse = ",")) %>%
  dplyr::mutate(compound = stringr::str_remove_all(compound, "[c()\\\"]")) %>%
  tidyr::pivot_wider(names_from = reactionSide,
                     values_from = compound) %>%
  dplyr::inner_join(rhea_table %>%
                      dplyr::select(rheaid, reactionEquation),
                    .,
                    by = "rheaid") %>%
  dplyr::distinct_all()


# Write table to output
readr::write_tsv(rhea_table, paste0("data/", run_date, "_rhea_db_parsed.tsv"))
readr::write_tsv(rhea_reactions, paste0("data/", run_date, "_rhea_db_reactions.tsv"))
