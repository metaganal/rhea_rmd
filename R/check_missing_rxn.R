library(argparser)
library(tidyverse)


# Define argument parser
p <- arg_parser('Check Missing RXN')
p <- add_argument(p, 'rxn_list',
                  help = 'List of available RXN files')
p <- add_argument(p, 'mol_list',
                  help = 'List of available MOL files')
p <- add_argument(p, 'reaction_table',
                  help = 'TSV file of RHEA reactions')
p <- add_argument(p, 'generic_table',
                  help = 'TSV file of RHEA generics reactive part CHEBI')
p <- add_argument(p, 'rxn_dir',
                  help = 'Directory of RXN files')
p <- add_argument(p, 'mol_dir',
                  help = 'Directory of MOL files')
p <- add_argument(p, 'out_dir',
                  help = 'Output directory of list and tables')

# Read arguments
ar <- parse_args(p)


# Check RXN list against RHEA reactions
rxn_list <- read_lines(ar$rxn_list) %>%
  gsub(pattern = '.rxn', replacement = '', fixed = T) %>%
  as.integer()
reactions <- read_tsv(ar$reaction_table) %>%
  mutate(L = gsub(':', '_', L, fixed = T),
         R = gsub(':', '_', R, fixed = T))
rhea_ids <- gsub('RHEA:', '', reactions$rheaid, fixed = T) %>% as.integer()
left_rxn <- rxn_list[rxn_list %in% (rhea_ids + 1)]
right_rxn <- sort(c(rxn_list[rxn_list %in% (rhea_ids + 2)], 28622))
no_rxn <- paste0('RHEA:', setdiff(rhea_ids, (left_rxn - 1)))


# List of reactions with no RXN file
no_rxn_reactions <- reactions %>%
  filter(rheaid %in% no_rxn)
compound_to_check <- c('CHEBI_13193',
                       'CHEBI_17499') %>%
  paste(collapse = '|')
compound_to_filterout <- c('CHEBI_13193[0-9]',
                           'CHEBI_17499[0-9]') %>%
  paste(collapse = '|')
a_reactions <- reactions %>%
  filter(grepl(compound_to_check, L)|
           grepl(compound_to_check, R)) %>%
  filter(!grepl(compound_to_filterout, L)|
           !grepl(compound_to_filterout, R))


# Check list of molfiles
molfiles <- read_lines(ar$mol_list)


# Prepare GENERIC to CHEBI conversion table
generic_raw <- read_tsv(ar$generic_table) %>%
  rename(generic = compoundAc) %>%
  mutate(generic = gsub(':', '_', generic, fixed = T),
         chebi = gsub(':', '_', chebi, fixed = T))
generic_mult_chebi <- generic_raw %>%
  count(generic) %>%
  filter(n > 1)
generic_chebi <- generic_raw %>%
  group_by(generic) %>%
  summarise(chebi = paste(chebi, collapse = ','))


# Prepare equation with no A or AH2
no_generic_rxn <- a_reactions %>%
  filter(!grepl('GENERIC', L), !grepl('GENERIC', R))
removed_a_reactions <- no_generic_rxn %>%
  mutate(L = str_split(L, ', '),
         R = str_split(R, ', ')) %>%
  unnest(L) %>%
  unnest(R) %>%
  filter(!grepl('CHEBI_13193$|CHEBI_17499$', L),
         !grepl('CHEBI_13193$|CHEBI_17499$', R)) %>%
  group_by(rheaid, reactionEquation) %>%
  summarise(n_l = length(unique(L)),
            n_r = length(unique(R)),
            L = paste(unique(L), collapse = ','),
            R = paste(unique(R), collapse = ',')) %>%
  ungroup()


# Convert GENERIC to CHEBI
# Extract reactions with GENERIC
generic_rxn <- a_reactions %>%
  filter(grepl('GENERIC', L) | grepl('GENERIC', R))


# Extract GENERIC converted to multiple CHEBI
generic_rxn_mult_chebi <- generic_rxn %>%
  filter(grepl(paste(generic_mult_chebi$generic, collapse = '|'), L) |
         grepl(paste(generic_mult_chebi$generic, collapse = '|'), R))


# Prepare equation with no A or AH2 and GENERIC converted to CHEBI
generic_rxn_no_a <- generic_rxn %>%
  mutate(L = str_split(L, ', '),
         R = str_split(R, ', ')) %>%
  unnest(L) %>%
  unnest(R) %>%
  filter(!grepl('CHEBI_13193$|CHEBI_17499$', L),
         !grepl('CHEBI_13193$|CHEBI_17499$', R)) %>%
  left_join(., generic_chebi, by = c('L' = 'generic')) %>%
  rename(chebi_l = chebi) %>%
  left_join(., generic_chebi, by = c('R' = 'generic')) %>%
  rename(chebi_r = chebi) %>%
  mutate(L = case_when(!is.na(chebi_l) ~ chebi_l,
                       TRUE ~ L),
         R = case_when(!is.na(chebi_r) ~ chebi_r,
                       TRUE ~ R)) %>%
  select(-chebi_l, -chebi_r) %>%
  group_by(rheaid, reactionEquation) %>%
  summarise(n_l = length(unique(L)),
            n_r = length(unique(R)),
            L = paste(unique(L), collapse = ','),
            R = paste(unique(R), collapse = ',')) %>%
  ungroup()


# Check if MOL files are available for RXN making
chebi_needed <- bind_rows(generic_rxn_no_a, removed_a_reactions) %>%
  mutate(compounds = paste(L, ',', R),
         compounds = gsub(' ', '', compounds),
         compounds = str_split(compounds, ',')) %>%
  unnest(compounds) %>%
  pull(compounds) %>%
  unique()
molfiles_needed <- paste0(chebi_needed, '.mol')
missing_molfiles <- setdiff(molfiles_needed, molfiles)
if (length(missing_molfiles) > 0) {
  message(paste0('There are missing MOL file(s), logged in ',
                 ar$out_dir, '/missing_molfile.log'))
  write_lines(missing_molfiles, paste0(ar$out_dir, 'missing_molfile.log'))
}


# Final reaction table for RXN file making
rxn_table <- bind_rows(generic_rxn_no_a, removed_a_reactions) %>%
  mutate(rheaid = gsub('RHEA:', '', rheaid),
         L = gsub(' ', '', L),
         R = gsub(' ', '', R)) %>%
  select(-reactionEquation)
write_tsv(rxn_table, paste0(ar$out_dir, '/rhea_make_rxn.tsv'))
rheaid_for_rxn <- bind_rows(generic_rxn_no_a, removed_a_reactions) %>%
  pull(rheaid)
diff_reactions <- setdiff(no_rxn_reactions$rhea_id, rheaid_for_rxn)
write_lines(diff_reactions, paste0(ar$out_dir, '/no_rxn_rheaid.list'))


ltor_params <- rxn_table %>%
  mutate(rheaid = as.integer(rheaid) + 1) %>%
  pmap_chr(function(rheaid, L, R, ...) {
    paste(rheaid,
          '-m', mol_dir,
          '-o', rxn_dir,
          '-r', L,
          '-p', R)
  })
rtol_params <- rxn_table %>%
  mutate(rheaid = as.integer(rheaid) + 2) %>%
  pmap_chr(function(rheaid, L, R, ...) {
    paste(rheaid,
          '-m', mol_dir,
          '-o', rxn_dir,
          '-r', R,
          '-p', L)
  })
write_lines(ltor_params,
            paste0(ar$out_dir, '/left_rxn.txt'))
write_lines(rtol_params,
            paste0(ar$out_dir, '/right_rxn.txt'))
write_lines(c(ltor_params, rtol_params),
            paste0(ar$out_dir, '/rxn_params.txt'))
