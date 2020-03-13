setwd('~/rhea_rmd')
set.seed(0)


make_network <- function(db, net_seed, n_step = 1, current_step = 0) {
  if (is.character(net_seed) && length(net_seed) == 1) {
    network <- db %>%
      dplyr::filter(compoundID_L == net_seed |
                      compoundID_R == net_seed) %>%
      dplyr::distinct_all()
    current_step <- current_step + 1
  } else if (is.character(net_seed) && length(net_seed) > 1) {
    network <- db %>%
      dplyr::filter(compoundID_L %in% net_seed |
                      compoundID_R %in% net_seed) %>%
      dplyr::distinct_all()
    current_step <- current_step + 1
  }
  nodes <- tibble(id = c(network$compoundID_L,
                         network$compoundID_R),
                  label = "",
                  mask = c(network$compoundName_L,
                            network$compoundName_R)) %>%
    dplyr::distinct_all() %>%
    dplyr::arrange(id)
  if (current_step == n_step) {
    directions <- db %>%
      dplyr::filter(rheaid %in% network$rheaid) %>%
      dplyr::select(rheaid, direction)
    edges_tmp <- network %>%
      dplyr::select(rheaid, compoundID_L, compoundID_R) %>%
      dplyr::inner_join(., directions, by = 'rheaid')
    edges_reverse <- edges_tmp %>%
      dplyr::filter(direction == 'L<>R') %>%
      dplyr::rename(from = compoundID_R,
                    to = compoundID_L,
                    label = rheaid) %>%
      dplyr::select(from, to, label, direction)
    edges <- edges_tmp %>%
      dplyr::rename(from = compoundID_L,
                    to = compoundID_R,
                    label = rheaid) %>%
      dplyr::select(from, to, label, direction) %>%
      dplyr::bind_rows(., edges_reverse)
    return(list(edge = edges,
                node = nodes))
  } else {
    make_network(db = db,
                 net_seed = nodes %>% pull(id),
                 n_step = n_step,
                 current_step = current_step)
  }
}


library(tidyverse)
library(visNetwork)
library(tidygraph)
library(ggraph)


ar <- commandArgs(trailingOnly = TRUE)


# Usage
if (length(ar) < 4) {
  usage <- paste('Usage: Rscript generate_rhea_graph.R',
                 '[seed_CHEBI] [network_step] [out_graphml] [out_pdf]')
  stop(usage)
}


rhea_reactions <- read_tsv('data/rhea_direction_annotated_reactions.tsv')


rhea_network <- make_network(db = rhea_reactions,
                             net_seed = ar[1],
                             n_step = ar[2])


rhea_graph <- tbl_graph(nodes = rhea_network$node,
                        edges = rhea_network$edge,
                        directed = TRUE)


visNetwork(nodes = rhea_network$node,
           edges = rhea_network$edge) %>%
  visNodes() %>%
  visEdges(arrows = "to") %>%
  visInteraction(hover = T) %>%
  visEvents(selectNode  = "function(e){
            var label_info = this.body.data.nodes.get({
            fields: ['mask', 'label'],
            filter: function (item) {
            return item.id === e.node
            },
            returnType :'Array'
            });
            this.body.data.nodes.update({id: e.node, mask : label_info[0].label, label : label_info[0].mask});
            }") %>%
  visEvents(blurNode  = "function(e){
            var label_info = this.body.data.nodes.get({
            fields: ['mask', 'label'],
            filter: function (item) {
            return item.id === e.node
            },
            returnType :'Array'
            });
            this.body.data.nodes.update({id: e.node, mask : label_info[0].label, label : label_info[0].mask});
  }")



igraph::write_graph(graph = rhea_graph,
                    file = ar[3],
                    format = 'graphml')


graph_vis <- ggraph(rhea_graph) +
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')),
                 end_cap = circle(3, 'mm')) +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = mask), repel = TRUE) +
  theme_graph()

ggplot2::ggsave(filename = ar[4],
                plot = graph_vis, width = 8, height = 8)

