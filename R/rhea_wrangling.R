# Set working directory to the git repository
setwd("~/rhea_rmd")


# Load package
library(tidyverse)


# Read RHEA table
dframe <- read_tsv("~/digzyme/rhea_20200313.tsv")


# Extract reactions
reaction_table <- dframe %>%
  mutate(reactionSide = gsub(".*[0-9]{5}_", "", reactionSide)) %>%
  select(rheaid, reactionEquation, reactionSide, compoundID, compoundName) %>%
  arrange(rheaid, reactionSide)


# Extract reaction compounds
rhea_table <- reaction_table %>%
  group_by(rheaid, reactionSide) %>%
  summarise(compound = paste(list(compoundID), collapse = ",")) %>%
  mutate(compound = str_remove_all(compound, "[c()\\\"]")) %>%
  pivot_wider(names_from = reactionSide, values_from = compound) %>%
  inner_join(reaction_table %>% select(rheaid, reactionEquation), .,
             by = "rheaid") %>%
  distinct_all()


# Write table to output
write_tsv(rhea_table, "data/rhea_rheaid_participants.tsv")
write_tsv(reaction_table, "data/rhea_id_side_compound.tsv")
