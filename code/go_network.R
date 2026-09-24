# =============================================================================
# E-GO-D network: exposure-GO-disease network in Cytoscape
# =============================================================================

rm(list = ls())

library(RCy3)
library(tidyverse)
library(writexl)

# -----------------------------------------------------------------------------
# Input (load before running):
#   EGoD_edges.csv: edge table exported from the ExposomeX platform
#   (columns: Source, Target, "Source group", "Target group",
#    "Source preferred name", "Target preferred name", Database)
# -----------------------------------------------------------------------------
# edges_raw <- read_csv("EGoD_edges.csv")

# Pollutant abbreviation <-> ExposomeX exposure ID
pollutant_mapping <- data.frame(
  source_id = c("EX:E000005038", "EX:E000060919", "EX:E000005055", "EX:E000005278",
                "EX:E000005052", "EX:E000004290", "EX:E000065012", "EX:E000002599",
                "EX:E000004331", "EX:E028250140", "EX:E028250143", "EX:E028250141",
                "EX:E028245645", "EX:E051947117", "EX:E002558306", "EX:E000191299",
                "EX:E000098621", "EX:E051947118", "EX:E000105427", "EX:E023920745",
                "EX:E024269009", "EX:E090574179", "EX:E000166282", "EX:E008359588",
                "EX:E008578608", "EX:E083973165", "EX:E009759597", "EX:E000003163",
                "EX:E000078878", "EX:E000023751", "EX:E000020455", "EX:E000002567",
                "EX:E000002860", "EX:E000003780", "EX:E000026053", "EX:E000002109",
                "EX:E090564678", "EX:E000011447", "EX:E000006270", "EX:E000006264",
                "EX:E000538902", "EX:E000558126", "EX:E000008218", "EX:E000006630",
                "EX:E000088028", "EX:E000068871", "EX:E000062971", "EX:E000090236",
                "EX:E000075021", "EX:E000020066", "EX:E000018334", "EX:E000142364",
                "EX:E013156095", "EX:E000015745", "EX:E000019493", "EX:E011214727",
                "EX:E000021564", "EX:E000022460", "EX:E000022401", "EX:E000096046",
                "EX:E000022448", "EX:E000022464", "EX:E004815849", "EX:E005270968",
                "EX:E004815779", "EX:E000022407", "EX:E000022444", "EX:E005271068",
                "EX:E004815817", "EX:E004811303"),
  abbreviation = c("SFM", "4AN6C2B", "SFSX", "TMTH", "SFTHI", "NOR", "ENOR", "CIP",
                   "OFLX", "TET", "OTC", "AM", "DOXY", "IMI", "NIT", "ACE", "THM",
                   "CLO", "THI", "IPP", "CYC", "DIN", "IMID", "DM-THM", "DN-IMI",
                   "5-OH-IMI", "DM-ACE", "FP", "DMM", "CBDZ", "CPS-O", "DZN", "DDVP",
                   "MALT", "3OH-CBF", "ATZ", "ATZ-DE-DSI-2HD", "BPF", "BPS", "BPC",
                   "BPE", "BPP", "2-OHN", "1-OHN", "3-OHF", "2-OHF", "2-OHP", "1-OHP",
                   "4-OHP", "1-OHPYR", "3-PBA", "4-F-3-PBA", "2-MB-3CA", "trans-DCCA",
                   "cis-DCCA", "IMI-olefin", "3,5,6-TCP", "V", "Fe", "Co", "Cu", "Zn",
                   "As", "Se", "Sr", "Mo", "Cd", "Te", "Tl", "Pb")
)

replace_with_abbreviation <- function(source_id) {
  if (is.na(source_id)) return(NA)
  idx <- match(source_id, pollutant_mapping$source_id)
  if (!is.na(idx)) {
    return(pollutant_mapping$abbreviation[idx])
  } else {
    return(source_id)
  }
}

# -----------------------------------------------------------------------------
# Clean edge data
# -----------------------------------------------------------------------------
print(names(edges_raw))

# Replace source preferred name with abbreviation
edges_raw <- edges_raw %>%
  mutate(
    `Source preferred name` = sapply(Source, replace_with_abbreviation)
  )

# Unify disease nodes to Hypertension
edges_raw <- edges_raw %>%
  mutate(
    Target = as.character(Target),
    Target = ifelse(`Target group` == "Disease", "Hypertension", Target),
    `Target preferred name` = ifelse(`Target group` == "Disease", "Hypertension", `Target preferred name`)
  )

# -----------------------------------------------------------------------------
# Build node table
# -----------------------------------------------------------------------------
nodes_source <- edges_raw %>%
  select(id = Source,
         label = `Source preferred name`,
         group = `Source group`) %>%
  distinct(id, .keep_all = TRUE)

nodes_target <- edges_raw %>%
  select(id = Target,
         label = `Target preferred name`,
         group = `Target group`) %>%
  distinct(id, .keep_all = TRUE)

Nodes_GO <- bind_rows(nodes_source, nodes_target) %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(
    id = ifelse(group == "Disease", "Hypertension", id),
    label = ifelse(group == "Disease", "Hypertension", label),
    color = case_when(
      group == "Exposure" ~ "chemical",
      group == "GO"       ~ "go",
      group == "Disease"  ~ "disease",
      TRUE                ~ "other"
    )
  ) %>%
  distinct(id, .keep_all = TRUE)

# -----------------------------------------------------------------------------
# Build edge table
# -----------------------------------------------------------------------------
Edges_GO <- edges_raw %>%
  mutate(
    source      = Source,
    target      = Target,
    interaction = "association",
    source.class = `Source group`,
    target.class = `Target group`,
    database    = Database,
    edge.color  = "#AAAAAA",
    type = case_when(
      target.class == "GO" ~ 1,
      target.class == "Disease" ~ 2,
      TRUE                      ~ 0
    )
  ) %>%
  select(source, target, interaction, source.class, target.class,
         database, edge.color, type) %>%
  distinct(source, target, .keep_all = TRUE)

print(paste("NA values in edge table:", sum(is.na(Edges_GO))))

# -----------------------------------------------------------------------------
# Node degree and size
# -----------------------------------------------------------------------------
degree_source <- Edges_GO %>% count(source, name = "n_source")
degree_target <- Edges_GO %>% count(target, name = "n_target")

node_degree <- full_join(degree_source, degree_target,
                         by = c("source" = "target")) %>%
  mutate(n = coalesce(n_source, 0) + coalesce(n_target, 0)) %>%
  select(id = source, n)

Nodes_GO <- Nodes_GO %>%
  left_join(node_degree, by = "id") %>%
  mutate(n = ifelse(is.na(n), 0, n))

Nodes_GO <- Nodes_GO %>%
  mutate(count = case_when(
    group == "Disease" ~ 140,
    group == "Exposure" ~ case_when(
      n >= 1 & n <= 5      ~ 10,
      n > 5  & n <= 15    ~ 30,
      n > 15 & n <= 30    ~ 50,
      n > 30 & n <= 50    ~ 70,
      n > 50 & n <= 90    ~ 90,
      n > 90 & n <= 200   ~ 110,
      n > 200             ~ 130,
      TRUE                ~ 5
    ),
    group == "GO" ~ case_when(
      n >= 1 & n <= 5      ~ 30,
      n > 5  & n <= 15    ~ 50,
      n > 15 & n <= 30    ~ 70,
      n > 30 & n <= 50    ~ 90,
      n > 50 & n <= 90    ~ 110,
      n > 90 & n <= 200   ~ 130,
      n > 200             ~ 150,
      TRUE                ~ 5
    ),
    TRUE ~ 20
  ))

print(paste("Nodes:", nrow(Nodes_GO)))
print(paste("Edges:", nrow(Edges_GO)))
print("Node groups:")
print(table(Nodes_GO$group))

# -----------------------------------------------------------------------------
# Create network in Cytoscape
# -----------------------------------------------------------------------------
cytoscapePing()
closeSession(F)

createNetworkFromDataFrames(
  nodes = Nodes_GO,
  edges = Edges_GO,
  title = "Exposure-GO-Hypertension Network",
  collection = "MyNetwork"
)

# -----------------------------------------------------------------------------
# Visual style
# -----------------------------------------------------------------------------
style_name <- "Style1"
defaults <- list(
  NODE_SHAPE = "ELLIPSE",
  NODE_SIZE = 20,
  EDGE_TRANSPARENCY = 200,
  EDGE_WIDTH = 1,
  EDGE_COLOR = "#AAAAAA",
  EDGE_TARGET_ARROW_SHAPE = "NONE",
  NODE_LABEL_POSITION = "S,W,c,0.00,0.00"
)

node_label <- mapVisualProperty('node label', 'label', 'p')
node_shape <- mapVisualProperty('Node Shape', 'group', 'd',
                                c("Exposure", "Disease", "GO"),
                                c("ELLIPSE", "ELLIPSE", "ELLIPSE"))
node_size <- mapVisualProperty('Node Size', 'count', 'd',
                               c(8, 20, 30, 50, 70, 90, 110, 120, 130, 140),
                               c(8, 20, 30, 50, 70, 90, 110, 120, 130, 140))
node_label_font_size <- mapVisualProperty('Node Label Font Size', 'group', 'd',
                                          c("Exposure", "Disease", "GO"),
                                          c(30, 30, 25))
node_z_location <- mapVisualProperty('Node Z Location', 'group', 'd',
                                     c("Exposure", "Disease", "GO"),
                                     c(3, 2, 1))

createVisualStyle(style_name,
                  defaults,
                  list(node_label,
                       node_shape,
                       node_size,
                       node_label_font_size,
                       node_z_location))

setVisualStyle(style_name)

setNodeColorMapping(
  mapping.type = 'd',
  table.column = 'group',
  table.column.values = c("Exposure", "Disease", "GO"),
  colors = c("#E74C3C", "#3498DB", "#A9A7D6"),
  style.name = style_name
)
