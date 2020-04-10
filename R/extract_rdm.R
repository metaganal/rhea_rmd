setwd('~/rhea_rmd')
set.seed(0)


library(tidyverse)


react_file <- 'data/rhea_direction_annotated_reactions.tsv'
rhea_reactions <- readr::read_tsv(react_file)
curated_reactions <- rhea_reactions %>%
  dplyr::select(rheaid, compoundID_L, compoundID_R, direction) %>%
  dplyr::slice(-grep('GENERIC', compoundID_L)) %>%
  dplyr::slice(-grep('GENERIC', compoundID_R)) %>%
  dplyr::slice(-grep(',', compoundID_L)) %>%
  dplyr::slice(-grep(',', compoundID_R)) %>%
  dplyr::mutate(compoundID_L = gsub(':', '_', compoundID_L),
                compoundID_R = gsub(':', '_', compoundID_R))

rhea_compounds <- dplyr::union(curated_reactions$compoundID_L,
                               curated_reactions$compoundID_R)


mol_files <- gsub('.mol', '', list.files(molfile_dir, './mol'))
usable_mol <- dplyr::intersect(mol_files, rhea_compounds)


curated_mol <- curated_reactions %>%
  dplyr::filter(compoundID_L %in% usable_mol,
                compoundID_R %in% usable_mol)


# Replace R# with R in mol files
# system(command = 'sed -i "s/R#/R/g" ~/digzyme/mol/*.mol')


# Convert MOLFILE to KCF
# dir.create('~/digzyme/kcf/')
# purrr::walk(.x = usable_mol,
#             .f = function(from) {
#               mol <- paste0('~/digzyme/mol/', from, '.mol')
#               kcf <- paste0('~/digzyme/kcf/', from, '.kcf')
#               shell <- paste('~/digzyme/kegg_tools/kcfco-1.1.3/makekcf',
#                              mol, kcf)
#               system(command = shell)
#             })


# dnode_atom <- tibble::tibble(d = 0,
#                              a = c(0, 1, 2)) %>%
#   dplyr::bind_rows(.,
#                    tibble::tibble(d = 1,
#                                   a = c(0, 1, 2)))
#
#
# # SIMCOMP between compound pair
# sim_out <- '~/digzyme/simcomp_results'
# dir.create(sim_out)
# purrr::walk2(.x = curated_mol$compoundID_L,
#              .y = curated_mol$compoundID_R,
#              .f = function(from, to) {
#                substrate <- paste0('~/digzyme/kcf/',
#                                    from, '.kcf')
#                product <- paste0('~/digzyme/kcf/',
#                                  to, '.kcf')
#                tools <- '~/digzyme/kegg_tools/simcomp-2.1.5/simcompp'
#                dir_out <- paste0(sim_out, '/', from, '_', to)
#                dir.create(dir_out)
#                purrr::walk2(.x = dnode_atom$d,
#                             .y = dnode_atom$a,
#                             .f = function(d, a,
#                                           f = substrate,
#                                           t = product) {
#                               kcf_out <- paste0(dir_out,
#                                                 '/d', d, '_a', a,
#                                                 '.kcf')
#                               shell <- paste(tools,
#                                              '--dnode', d, '--atom', a,
#                                              '--trick 0 --trickatom 0',
#                                              '--kcfoutput --output',
#                                              kcf_out, f, t)
#                               system(command = shell)
#                             })
#              })


# Extract RDM patterns from SIMCOMP results
# pl <- '~/digzyme/kegg_tools/rdmchk_m1.pl'
# bin <- '~/digzyme/kegg_tools/rdmchk-0.9.8.3.p/src/rdmchk'
# rdm_out <- '~/digzyme/rdmchk_results'
# dir.create(rdm_out)
# purrr::walk2(.x = curated_mol$compoundID_L,
#              .y = curated_mol$compoundID_R,
#              .f = function(from, to) {
#                substrate <- paste0('~/digzyme/kcf/',
#                                    from, '.kcf')
#                product <- paste0('~/digzyme/kcf/',
#                                  to, '.kcf')
#                simcomp_dir <- paste0(sim_out, '/', from, '_', to)
#                dir_out <- paste0(rdm_out, '/', from, '_', to)
#                dir.create(dir_out)
#                purrr::walk2(.x = dnode_atom$d,
#                             .y = dnode_atom$a,
#                             .f = function(d, a,
#                                           f = substrate,
#                                           t = product) {
#                               kcf_in <- paste0(simcomp_dir, '/d', d, '_a', a,
#                                               '.kcf')
#                               kcf_out <- paste0(dir_out, '/d', d, '_a', a,
#                                                 '_rdm.kcf')
#                               shell <- paste('perl', pl, bin,
#                                              kcf_in, kcf_out)
#                               system(command = shell)
#                             })
#              })


# rdm_files <- list.files(path = rdm_out, full.names = TRUE, recursive = TRUE)
# readr::write_lines(x = rdm_files, path = '~/digzyme/rdm_files.txt')
rdm_files <- read_lines('~/digzyme/rdm_files.txt')


parse_rdm <- function(rdm_file) {
  rdm_f <- read_lines(rdm_file)
  reaction <- basename(dirname(rdm_file))
  simcomp_param <- str_extract(rdm_file, 'd[01]_a[0-2]')

  if (length(rdm_f) == 0) {
    return(tibble::tibble(reaction = reaction,
                          simcomp_parameters = simcomp_param,
                          rdm = ''))
  }

  trash <- grep('ALIGN', rdm_f) - 1
  rdm_f <- rdm_f[1:trash]
  rdms <- str_remove(rdm_f[6:trash], "^[ 0-9]{0,100}")
  tibble::tibble(reaction = reaction,
                 simcomp_parameters = simcomp_param,
                 rdm = rdms)
}


rdm_dataframe <- rdm_files %>%
  purrr::map_dfr(parse_rdm) %>%
  dplyr::distinct(reaction, rdm, .keep_all = TRUE)
# failed rdmchk due to input file error (simcomp kcf) CHEBI_139218_CHEBI_17087
# https://www.rhea-db.org/reaction?id=26490 example of failed simcomp

# readr::write_tsv(rdm_dataframe, 'rhea_extracted_rdm.tsv')
# old_rdms <- readr::read_tsv('rhea_extracted_rdm.tsv')


# rdm_dataframe %>%
#   distinct(reaction, rdm) %>%
#   nrow()


multi_rdms <- rdm_dataframe %>%
  dplyr::group_by(reaction) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n > 1) %>%
  dplyr::pull(reaction)


multi_rdm_dataframe <- rdm_dataframe %>%
  dplyr::filter(reaction %in% multi_rdms)

single_rdm_dataframe <- rdm_dataframe %>%
  dplyr::filter(!reaction %in% multi_rdms)


curated_mol <- curated_mol %>%
  mutate(reaction = paste0(compoundID_L, '_', compoundID_R))


single_rdm_dataframe2 <- inner_join(single_rdm_dataframe,
                                    curated_mol %>%
                                      select(rheaid, reaction),
                                    by = 'reaction')

multi_rdm_dataframe2 <- inner_join(multi_rdm_dataframe,
                                    curated_mol %>%
                                      select(rheaid, reaction),
                                    by = 'reaction')


multi_rhea <- single_rdm_dataframe2 %>%
  group_by(reaction) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>%
  pull(reaction)
multi_rhea2 <- single_rdm_dataframe2 %>%
  filter(reaction %in% multi_rhea)


missing_rdm <- multi_rdm_dataframe2 %>%
  filter(rdm == "")


x <- single_rdm_dataframe2 %>%
  filter(rdm != "") %>%
  group_by(rdm) %>%
  nest() %>%
  glimpse()
  summarise(reactions = list(reaction))


y <- multi_rdm_dataframe2 %>%
  filter(rdm != "") %>%
  group_by(reaction) %>%
  nest(data = c(simcomp_parameters, rdm)) %>%
  glimpse()

y %>%
  group_by(data) %>%
  nest()

# Example of failed RDM https://www.rhea-db.org/reaction?id=21032,
# https://www.rhea-db.org/reaction?id=13577


save.image('extract_rdm.RData')
