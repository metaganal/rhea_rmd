setwd("~/rhea_rmd")

# Function to curate direction by directing compound, e.g. O2
curate_direction <- function(dframe, id) {
  substrate <- dframe %>%
    dplyr::slice(grep(id, L)) %>%
    dplyr::mutate(direction = "L>R")
  product <- dframe %>%
    dplyr::slice(grep(id, R)) %>%
    dplyr::mutate(direction = "L>R")
  directed <- dplyr::bind_rows(substrate, product)
}


# Load libraries
library(tidyverse)

run_date <- gsub(Sys.Date(), pattern = "-", replacement = "", fixed = T)


# Parse command line arguments
ar <- commandArgs(trailingOnly = TRUE)


# Read input TSV
rhea_reactions <- readr::read_tsv(ar[[1]])
rhea_table <- readr::read_tsv(ar[[2]])


# Calculate compound usage frequency
compound_freq <- rhea_table %>%
  dplyr::count(compoundID, compoundName) %>%
  dplyr::arrange(desc(n))


# Make vectors of CHEBI ID to curate with
chebi_to_curate <- c("CHEBI:15379", "CHEBI:30616", "CHEBI:456216",
                     "CHEBI:456215", "CHEBI:16526", "CHEBI:59789") %>%
  purrr::set_names(c("O2", "ATP", "ADP", "AMP", "CO2", "AMS"))


# Curate based on compounds presence
curated_dfs <- chebi_to_curate %>%
  purrr::map(function(chebi) {
    curate_direction(dframe = rhea_reactions, id = chebi)
  })


# Aggregate curated dataframes
directed_rhea <- purrr::reduce(.x = curated_dfs,
                               .f = dplyr::bind_rows) %>%
  dplyr::distinct(rheaid, .keep_all = TRUE) %>%
  dplyr::arrange(rheaid)


# Bidirectional reactions table
bidirectional_rhea <- rhea_reactions %>%
  dplyr::filter(!rheaid %in% directed_rhea$rheaid) %>%
  dplyr::mutate(direction = 'L<>R')


# Extract compounds with frequent usage
frequent_compound <- compound_freq %>%
  dplyr::filter(n > 10) %>%
  dplyr::pull(compoundID)



to_check <- rhea_table %>%
  dplyr::filter(rheaid %in% directed_rhea$rheaid,
                !compoundID %in% c("CHEBI:456215",
                                   "CHEBI:456216",
                                   "CHEBI:59789",
                                   "CHEBI:16526",
                                   "CHEBI:30616",
                                   "CHEBI:15379")) %>%
  dplyr::select(rheaid, reactionEquation) %>%
  dplyr::inner_join(., directed_rhea %>% dplyr::select(rheaid, direction),
                    by = 'rheaid') %>%
  dplyr::distinct(rheaid, .keep_all = TRUE)


# Extract reactions with no frequent compounds
main_directed <- rhea_table %>%
  dplyr::filter(rheaid %in% directed_rhea$rheaid,
                !compoundID %in% frequent_compound) %>%
  dplyr::group_by(rheaid, reactionSide, reactionEquation) %>%
  dplyr::summarise(compoundID = paste(compoundID, collapse = ','),
                   compoundName = paste(compoundName, collapse = ',')) %>%
  tidyr::pivot_wider(names_from = reactionSide,
                     values_from = c(compoundID, compoundName))


# Extract bidirectional reactions with no frequent compounds
main_bidirectional <- rhea_table %>%
  dplyr::filter(rheaid %in% bidirectional_rhea$rheaid,
                !compoundID %in% frequent_compound) %>%
  dplyr::group_by(rheaid, reactionSide, reactionEquation) %>%
  dplyr::summarise(compoundID = paste(compoundID, collapse = ','),
                   compoundName = paste(compoundName, collapse = ',')) %>%
  tidyr::pivot_wider(names_from = reactionSide,
                     values_from = c(compoundID, compoundName))


# Extract reactions with cofactors as main participants
cofactor_directed <- rhea_table %>%
  dplyr::filter(rheaid %in% directed_rhea$rheaid,
                compoundID %in% frequent_compound) %>%
  dplyr::group_by(rheaid, reactionSide, reactionEquation) %>%
  dplyr::summarise(compoundID = paste(compoundID, collapse = ','),
                   compoundName = paste(compoundName, collapse = ',')) %>%
  tidyr::pivot_wider(names_from = reactionSide,
                     values_from = c(compoundID, compoundName))


# Extract bidirectional reactions with cofactors as main participants
cofactor_bidirectional <- rhea_table %>%
  dplyr::filter(rheaid %in% bidirectional_rhea$rheaid,
                compoundID %in% frequent_compound) %>%
  dplyr::group_by(rheaid, reactionSide, reactionEquation) %>%
  dplyr::summarise(compoundID = paste(compoundID, collapse = ','),
                   compoundName = paste(compoundName, collapse = ',')) %>%
  tidyr::pivot_wider(names_from = reactionSide,
                     values_from = c(compoundID, compoundName))


# Extract reactions which involves only frequently used compounds
generic_main_directed <- cofactor_directed %>%
  dplyr::filter(!rheaid %in% main_directed$rheaid)
generic_main_bidirectional <- cofactor_bidirectional %>%
  dplyr::filter(!rheaid %in% main_bidirectional$rheaid)


#
cofactor_directed <- cofactor_directed %>%
  dplyr::filter(!rheaid %in% generic_main_directed$rheaid)
cofactor_bidirectional <- cofactor_bidirectional %>%
  dplyr::filter(!rheaid %in% generic_main_bidirectional$rheaid)


# Final curated RHEA reactions (mono for unidirectional reaction)
mono_reaction <- main_directed %>%
  tidyr::drop_na() %>%
  dplyr::slice(c(grep(',', compoundID_L, invert = TRUE)),
               grep(',', compoundID_R, invert = TRUE)) %>%
  dplyr::inner_join(.,
                    directed_rhea %>%
                      dplyr::select(rheaid, direction),
                    by = 'rheaid') %>%
  dplyr::distinct_all()
bi_reaction <- main_bidirectional %>%
  tidyr::drop_na() %>%
  dplyr::slice(c(grep(',', compoundID_L, invert = TRUE)),
               grep(',', compoundID_R, invert = TRUE)) %>%
  dplyr::filter(compoundID_L != compoundID_R) %>%
  dplyr::mutate(direction = 'L<>R') %>%
  dplyr::distinct_all()
full_reaction <- dplyr::bind_rows(mono_reaction, bi_reaction)


# Save curated dataframe to files
readr::write_tsv(main_directed, paste0("data/", run_date, "_rhea_reactants_unidirectional.tsv"))
readr::write_tsv(main_bidirectional, paste0("data/", run_date, "_rhea_reactants_bidirectional.tsv"))
readr::write_tsv(mono_reaction, paste0("data/", run_date, "_rhea_reactions_unidirectional.tsv"))
readr::write_tsv(bi_reaction, paste0("data/", run_date, "_rhea_reactions_bidirectional.tsv"))
readr::write_tsv(cofactor_directed, paste0("data/", run_date, "_rhea_cofactor_unidirectional.tsv"))
readr::write_tsv(cofactor_bidirectional, paste0("data/", run_date, "_rhea_cofactor_bidirectional.tsv"))
readr::write_tsv(generic_main_directed, paste0("data/", run_date, "_rhea_generic_unidirectional.tsv"))
readr::write_tsv(generic_main_bidirectional, paste0("data/", run_date, "_rhea_generic_bidirectional.tsv"))
readr::write_tsv(compound_freq, paste0("data/", run_date, "_rhea_compound_usage.tsv"))
readr::write_tsv(full_reaction, paste0("data/", run_date, "_rhea_reactions_annotated.tsv"))
