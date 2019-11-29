setwd("~/digzyme/")

# Function to curate direction by directing compound, e.g. O2
curate_direction <- function(dframe, id) {
  substrate <- dframe %>%
    dplyr::slice(grep(id, L)) %>%
    dplyr::mutate(direction = "L>R")
  product <- dframe %>%
    dplyr::slice(grep(id, R)) %>%
    dplyr::mutate(direction = "L<R")
  directed <- dplyr::bind_rows(substrate, product)
}


library(tidyverse)


dframe <- read_tsv("rhea_id_side_compound.tsv")


# Calculate compound usage frequency
compound_freq <- dframe %>%
  count(compoundID, compoundName) %>%
  arrange(desc(n))


rhea_table <- read_tsv("rhea_rheaid_participants.tsv")


# Make vectors of CHENI ID to curate with
chebi_to_curate <- c("CHEBI:15379", "CHEBI:30616", "CHEBI:456216",
                     "CHEBI:456215", "CHEBI:16526", "CHEBI:59789") %>%
  purrr::set_names(c("O2", "ATP", "ADP", "AMP", "CO2", "AMS"))


# Curate based on compounds presence
curated_dfs <- chebi_to_curate %>%
  purrr::map(function(chebi) {
    curate_direction(dframe = rhea_table, id = chebi)
  })


# Aggregate curated dataframes
directed_rhea <- reduce(.x = curated_dfs,
                        .f = dplyr::bind_rows) %>%
  dplyr::distinct(rheaid, .keep_all = TRUE) %>%
  dplyr::arrange(rheaid)


# Extract compounds with frequent usage
frequent_compound <- compound_freq %>%
  dplyr::filter(n > 10) %>%
  dplyr::pull(compoundID)



to_check <- dframe %>%
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
main_directed <- dframe %>%
  dplyr::filter(rheaid %in% directed_rhea$rheaid,
                !compoundID %in% frequent_compound) %>%
  dplyr::group_by(rheaid, reactionSide, reactionEquation) %>%
  dplyr::summarise(compoundID = paste(compoundID, collapse = ','),
                   compoundName = paste(compoundName, collapse = ',')) %>%
  tidyr::pivot_wider(names_from = reactionSide,
                     values_from = c(compoundID, compoundName))


# Extract reactions with cofactors as main participants
cofactor_directed <- dframe %>%
  dplyr::filter(rheaid %in% directed_rhea$rheaid,
                compoundID %in% frequent_compound) %>%
  dplyr::group_by(rheaid, reactionSide, reactionEquation) %>%
  dplyr::summarise(compoundID = paste(compoundID, collapse = ','),
                   compoundName = paste(compoundName, collapse = ',')) %>%
  tidyr::pivot_wider(names_from = reactionSide,
                     values_from = c(compoundID, compoundName))


generic_main_directed <- cofactor_directed %>%
  dplyr::filter(!rheaid %in% main_directed$rheaid)


#
cofactor_directed <- cofactor_directed %>%
  dplyr::filter(!rheaid %in% generic_main_directed$rheaid)


# Final curated RHEA reactions
mono_reaction <- main_directed %>%
  tidyr::drop_na() %>%
  dplyr::slice(c(grep(',', compoundID_L, invert = TRUE)),
               grep(',', compoundID_R, invert = TRUE)) %>%
  dplyr::distinct_all()


# Save curated dataframe to files
write_tsv(main_directed, "rhea_directed_reactants.tsv")
write_tsv(mono_reaction, "rhea_directed_mono_reactions.tsv")
write_tsv(cofactor_directed, "rhea_directed_cofactors.tsv")
write_tsv(generic_main_directed, "rhea_directed_generic.tsv")
write_tsv(compound_freq, "rhea_compound_usage.tsv")
