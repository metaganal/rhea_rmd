library(tidyverse)


ar <- commandArgs(trailingOnly = TRUE)
reaction_table <- ar[[1]]
generic_table <- ar[[2]]
compound_freq_chebi <- ar[[3]]
compound_freq_generic <- ar[[4]]


# Read reactions
reactions <- read_tsv(reaction_table) %>%
  mutate(L = gsub(':', '_', L, fixed = T),
         R = gsub(':', '_', R, fixed = T),
         L = gsub(' ', '', L, fixed = T),
         R = gsub(' ', '', R, fixed = T),
         L = str_split(L, ','),
         R = str_split(R, ',')) %>%
  unnest(L) %>%
  unnest(R)


# Prepare GENERIC to CHEBI conversion table
generic_raw <- read_tsv(generic_table) %>%
  rename(generic = compoundAc) %>%
  mutate(generic = gsub(':', '_', generic, fixed = T),
         chebi = gsub(':', '_', chebi, fixed = T))
generic_mult_chebi <- generic_raw %>%
  count(generic) %>%
  filter(n > 1)
generic_chebi <- generic_raw %>%
  group_by(generic) %>%
  summarise(chebi = paste(chebi, collapse = ','))


# Convert GENERIC to CHEBI
no_generic_reactions <- reactions %>%
  left_join(., generic_chebi, by = c('L' = 'generic')) %>%
  rename(chebi_l = chebi) %>%
  left_join(., generic_chebi, by = c('R' = 'generic')) %>%
  rename(chebi_r = chebi) %>%
  mutate(L = case_when(!is.na(chebi_l) ~ chebi_l,
                       TRUE ~ L),
         R = case_when(!is.na(chebi_r) ~ chebi_r,
                       TRUE ~ R)) %>%
  select(-chebi_l, -chebi_r)


# Calculate compound frequency of CHEBI converted reactions
substrate_frequency <- no_generic_reactions %>%
  select(rheaid, reactionEquation, L) %>%
  distinct_all() %>%
  count(L) %>%
  rename(compound = L, substrate_freq = n)
product_frequency <- no_generic_reactions %>%
  select(rheaid, reactionEquation, R) %>%
  distinct_all() %>%
  count(R) %>%
  rename(compound = R, product_freq = n)
compound_frequency <- full_join(substrate_frequency,
                                product_frequency,
                                by = 'compound') %>%
  modify_if(is.integer, function(x) if_else(is.na(x), 0L, x)) %>%
  arrange(compound)


# Calculate frequency with GENERIC as is
generic_substrate_frequency <- reactions %>%
  select(rheaid, reactionEquation, L) %>%
  distinct_all() %>%
  count(L) %>%
  rename(compound = L, substrate_freq = n)
generic_product_frequency <- reactions %>%
  select(rheaid, reactionEquation, R) %>%
  distinct_all() %>%
  count(R) %>%
  rename(compound = R, product_freq = n)
generic_compound_frequency <- full_join(substrate_frequency,
                                        product_frequency,
                                        by = 'compound') %>%
  modify_if(is.integer, function(x) if_else(is.na(x), 0L, x)) %>%
  arrange(compound)


# Output to file
dir.create(dirname(compound_freq_chebi), recursive = TRUE)
dir.create(dirname(compound_freq_generic), recursive = TRUE)
write_tsv(compound_frequency, compound_freq_chebi)
write_tsv(generic_compound_frequency, compound_freq_generic)
