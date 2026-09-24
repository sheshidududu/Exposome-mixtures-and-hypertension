# =============================================================================
# One-click network visualisation (SG core features only)
# Only the features retained in the final SG models (Is_Final == "Yes")
# Node fill colour by feature type; edges light grey
# =============================================================================

rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, igraph, RCy3, writexl, ggplot2, cowplot)

tryCatch({
  cytoscapePing()
}, error = function(e) {
  stop("Cannot connect to Cytoscape; make sure Cytoscape is running.")
})

# -----------------------------------------------------------------------------
# Input (load before running):
#   res_go      <- readRDS("GO/eSet_result.rds")       # E-GO-D results
#   res_gene    <- readRDS("gene153/eSet_result.rds")  # E-P-D results
#   res_pathway <- readRDS("pathway_v4/eSet_result.rds")  # E-Pa-D results
#   go_matrix      <- read.csv("pollutant_go_matrix.csv", row.names = 1,
#                              check.names = FALSE)
#   gene_matrix    <- read.csv("pollutant_gene_matrix.csv", row.names = 1,
#                              check.names = FALSE)
#   pathway_matrix <- read.csv("pollutant_pathway_matrix.csv", row.names = 1,
#                              check.names = FALSE)
# -----------------------------------------------------------------------------

imp_go      <- res_go$base_imp$Exposome
imp_gene    <- res_gene$base_imp$Exposome
imp_pathway <- res_pathway$base_imp$Exposome

cat("GO importance rows:", nrow(imp_go), "\n")
cat("Gene importance rows:", nrow(imp_gene), "\n")
cat("Pathway importance rows:", nrow(imp_pathway), "\n")

base_all <- bind_rows(
  imp_go      %>% select(Feature, Importance_pct),
  imp_gene    %>% select(Feature, Importance_pct),
  imp_pathway %>% select(Feature, Importance_pct)
) %>%
  mutate(Feature = paste0("Exposome_", Feature)) %>%
  group_by(Feature) %>%
  summarise(Base_Imp = max(Importance_pct, na.rm = TRUE), .groups = "drop")

cat("Merged base_all rows:", nrow(base_all), "\n")
if (nrow(base_all) == 0) {
  stop("base_all is empty; check base_imp$Exposome in the modelling results.")
}
print(head(base_all, 3))

# Features retained in the final SG models
final_go      <- res_go$sg$imp_df      %>% filter(VarType == "Exposome") %>% pull(Features)
final_gene    <- res_gene$sg$imp_df    %>% filter(VarType == "Exposome") %>% pull(Features)
final_pathway <- res_pathway$sg$imp_df %>% filter(VarType == "Exposome") %>% pull(Features)
final_feats   <- unique(c(final_go, final_gene, final_pathway))
cat("Final SG feature count:", length(final_feats), "\n")

# -----------------------------------------------------------------------------
# Pollutant groups and colors
# -----------------------------------------------------------------------------
groups_fine <- list(
  Sulfonamides     = c("SFM", "4AN6C2B", "SFSX", "TMTH", "SFTHI"),
  Fluoroquinolones = c("NOR", "ENOR", "CIP", "OFLX"),
  Tetracyclines    = c("TET", "OTC", "AM", "DOXY"),
  Neonicotinoids   = c("IMI", "NIT", "ACE", "THM", "CLO", "THI",
                       "IPP", "CYC", "DIN", "IMID", "DM-THM", "DN-IMI",
                       "5-OH-IMI", "DM-ACE", "IMI-olefin"),
  Organophosphates = c("DZN", "CPS-O", "3,5,6-TCP", "DDVP", "MALT"),
  Carbamates       = c("3OH-CBF"),
  Pyrethroids      = c("3-PBA", "4-F-3-PBA", "2-MB-3CA",
                       "trans-DCCA", "cis-DCCA"),
  Herbicides       = c("ATZ", "ATZ-DE-DSI-2HD", "CBDZ", "DMM"),
  Phenylpyrazoles  = c("FP"),
  Bisphenols       = c("BPF", "BPS", "BPC", "BPE", "BPP"),
  OH_PAHs          = c("2-OHN", "1-OHN", "3-OHF", "2-OHF",
                       "2-OHP", "1-OHP", "4-OHP", "1-OHPYR"),
  Metals           = c("V", "Fe", "Co", "Cu", "Zn", "As", "Se",
                       "Sr", "Mo", "Cd", "Te", "Tl", "Pb")
)

user_colors <- c(
  "Sulfonamides"     = "#C0392B",
  "Fluoroquinolones" = "#E74C3C",
  "Tetracyclines"    = "#E67E22",
  "Neonicotinoids"   = "#F39C12",
  "Organophosphates" = "#F1C40F",
  "Carbamates"       = "#D4AC0D",
  "Pyrethroids"      = "#82E0AA",
  "Herbicides"       = "#58D68D",
  "Phenylpyrazoles"  = "#2ECC71",
  "Bisphenols"       = "#5DADE2",
  "OH_PAHs"          = "#3498DB",
  "Metals"           = "#1A5276"
)

pollutant_map <- data.frame(
  RawName = unlist(groups_fine),
  Group   = rep(names(groups_fine), lengths(groups_fine)),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# Node table
# -----------------------------------------------------------------------------
nodes_cy <- base_all %>%
  mutate(
    id = Feature,
    Is_Final = ifelse(id %in% final_feats, "Yes", "No"),
    Plot_Size = sqrt(Base_Imp)
  ) %>%
  mutate(
    Type = case_when(
      str_starts(id, "Exposome_GO_")      ~ "GO",
      str_starts(id, "Exposome_Pathway_") ~ "Pathway",
      str_starts(id, "Exposome_Gene_")    ~ "Protein",
      TRUE                                ~ "Pollutant"
    ),
    Label = case_when(
      Type == "GO"      ~ str_remove(id, "^Exposome_GO_"),
      Type == "Pathway" ~ str_remove(id, "^Exposome_Pathway_"),
      Type == "Protein" ~ str_remove(id, "^Exposome_Gene_"),
      Type == "Pollutant" ~ str_remove(id, "^Exposome_"),
      TRUE ~ id
    ),
    RawName = str_remove(id, "^Exposome_")
  )

# Match pollutant class
nodes_cy <- nodes_cy %>%
  left_join(pollutant_map, by = "RawName") %>%
  mutate(Type = ifelse(Type == "Pollutant" & !is.na(Group), Group, Type)) %>%
  select(-RawName, -Group)

# Keep only SG final core features
nodes_cy <- nodes_cy %>% filter(Is_Final == "Yes")

# Add disease node
if (!"Hypertension" %in% nodes_cy$id) {
  nodes_cy <- bind_rows(
    nodes_cy,
    tibble(id = "Hypertension", Base_Imp = NA, Is_Final = "Yes",
           Plot_Size = if (nrow(nodes_cy) > 0) max(nodes_cy$Plot_Size, na.rm = TRUE) * 1.5 else 30,
           Type = "Disease", Label = "Hypertension")
  )
} else {
  max_size <- max(nodes_cy$Plot_Size, na.rm = TRUE)
  nodes_cy <- nodes_cy %>%
    mutate(Plot_Size = ifelse(id == "Hypertension", max_size * 1.5, Plot_Size))
}

cat("Nodes after filtering:", nrow(nodes_cy), "\n")
cat("Node counts by type:\n")
print(table(nodes_cy$Type))

# -----------------------------------------------------------------------------
# Standardize matrix row names to "EX:E" format
# -----------------------------------------------------------------------------
standardize_rownames <- function(mat) {
  rownames(mat) <- trimws(rownames(mat))
  rownames(mat) <- gsub("^X", "", rownames(mat))
  rownames(mat) <- gsub("^EX", "EX", rownames(mat))
  rownames(mat) <- ifelse(!startsWith(rownames(mat), "EX:"),
                          paste0("EX:", rownames(mat)),
                          rownames(mat))
  return(mat)
}

go_matrix      <- standardize_rownames(go_matrix)
gene_matrix    <- standardize_rownames(gene_matrix)
pathway_matrix <- standardize_rownames(pathway_matrix)

# -----------------------------------------------------------------------------
# Pollutant abbreviation -> EXID mapping
# -----------------------------------------------------------------------------
pollutant_mapping <- c(
  "SFM" = "EX:E000005038", "4AN6C2B" = "EX:E000060919", "SFSX" = "EX:E000005055",
  "TMTH" = "EX:E000005278", "SFTHI" = "EX:E000005052", "NOR" = "EX:E000004290",
  "ENOR" = "EX:E000065012", "CIP" = "EX:E000002599", "OFLX" = "EX:E000004331",
  "TET" = "EX:E028250140", "OTC" = "EX:E028250143", "AM" = "EX:E028250141",
  "DOXY" = "EX:E028245645", "IMI" = "EX:E051947117", "NIT" = "EX:E002558306",
  "ACE" = "EX:E000191299", "THM" = "EX:E000098621", "CLO" = "EX:E051947118",
  "THI" = "EX:E000105427", "IPP" = "EX:E023920745", "CYC" = "EX:E024269009",
  "DIN" = "EX:E090574179", "IMID" = "EX:E000166282", "DM-THM" = "EX:E008359588",
  "DN-IMI" = "EX:E008578608", "5-OH-IMI" = "EX:E083973165", "DM-ACE" = "EX:E009759597",
  "FP" = "EX:E000003163", "DMM" = "EX:E000078878", "CBDZ" = "EX:E000023751",
  "CPS-O" = "EX:E000020455", "DZN" = "EX:E000002567", "DDVP" = "EX:E000002860",
  "MALT" = "EX:E000003780", "3OH-CBF" = "EX:E000026053", "ATZ" = "EX:E000002109",
  "ATZ-DE-DSI-2HD" = "EX:E090564678", "BPF" = "EX:E000011447", "BPS" = "EX:E000006270",
  "BPC" = "EX:E000006264", "BPE" = "EX:E000538902", "BPP" = "EX:E000558126",
  "2-OHN" = "EX:E000008218", "1-OHN" = "EX:E000006630", "3-OHF" = "EX:E000088028",
  "2-OHF" = "EX:E000068871", "2-OHP" = "EX:E000062971", "1-OHP" = "EX:E000090236",
  "4-OHP" = "EX:E000075021", "1-OHPYR" = "EX:E000020066", "3-PBA" = "EX:E000018334",
  "4-F-3-PBA" = "EX:E000142364", "2-MB-3CA" = "EX:E013156095",
  "trans-DCCA" = "EX:E000015745", "cis-DCCA" = "EX:E000019493",
  "IMI-olefin" = "EX:E011214727", "3,5,6-TCP" = "EX:E000021564",
  "V" = "EX:E000022460", "Fe" = "EX:E000022401", "Co" = "EX:E000096046",
  "Cu" = "EX:E000022448", "Zn" = "EX:E000022464", "As" = "EX:E004815849",
  "Se" = "EX:E005270968", "Sr" = "EX:E004815779", "Mo" = "EX:E000022407",
  "Cd" = "EX:E000022444", "Te" = "EX:E005271068", "Tl" = "EX:E004815817",
  "Pb" = "EX:E004811303"
)

# Pollutant lookup table (core pollutants only)
pollutant_nodes <- nodes_cy %>%
  filter(Type %in% names(user_colors)) %>%
  mutate(Abbr = str_remove(id, "^Exposome_")) %>%
  select(id, Abbr)

pollutant_lookup <- pollutant_nodes %>%
  mutate(EXID = pollutant_mapping[Abbr]) %>%
  filter(!is.na(EXID))

cat("\nCore pollutant nodes:", nrow(pollutant_nodes), "\n")
cat("Pollutants with valid EXID:", nrow(pollutant_lookup), "\n")

# -----------------------------------------------------------------------------
# Edge generation from biological matrices
# -----------------------------------------------------------------------------
add_edges_from_matrix <- function(matrix_obj, target_prefix, edge_type) {
  if (target_prefix == "GO") {
    target_labels <- nodes_cy %>% filter(Type == "GO") %>% pull(Label)
    id_prefix <- "GO"
  } else if (target_prefix == "Pathway") {
    target_labels <- nodes_cy %>% filter(Type == "Pathway") %>% pull(Label)
    id_prefix <- "Pathway"
  } else if (target_prefix == "Protein") {
    target_labels <- nodes_cy %>% filter(Type == "Protein") %>% pull(Label)
    id_prefix <- "Gene"
  } else {
    return(tibble(source = character(), target = character(),
                  interaction = character()))
  }

  matrix_obj <- as.matrix(matrix_obj)
  mode(matrix_obj) <- "numeric"

  valid_rows <- intersect(rownames(matrix_obj), pollutant_lookup$EXID)
  valid_cols <- intersect(colnames(matrix_obj), target_labels)

  if (length(valid_rows) == 0 || length(valid_cols) == 0) {
    return(tibble(source = character(), target = character(),
                  interaction = character()))
  }

  edges <- tibble(source = character(), target = character(),
                  interaction = character())
  for (exid in valid_rows) {
    pol_node <- pollutant_lookup %>% filter(EXID == exid) %>% pull(id)
    if (length(pol_node) == 0) next
    for (col in valid_cols) {
      val <- matrix_obj[exid, col]
      if (!is.na(val) && val == 1) {
        target_id <- paste0("Exposome_", id_prefix, "_", col)
        edges <- bind_rows(edges, tibble(
          source = pol_node[1],
          target = target_id,
          interaction = edge_type
        ))
      }
    }
  }
  return(edges)
}

edges_go      <- add_edges_from_matrix(go_matrix,      "GO",      "activates")
edges_pathway <- add_edges_from_matrix(pathway_matrix, "Pathway", "activates")
edges_gene    <- add_edges_from_matrix(gene_matrix,    "Protein", "interacts_with")

edges_tmp <- bind_rows(edges_go, edges_pathway, edges_gene)
linked_pollutants <- if (nrow(edges_tmp) > 0) unique(edges_tmp$source) else character(0)

# Biological nodes -> disease
edges_bio <- nodes_cy %>%
  filter(Type %in% c("GO", "Pathway", "Protein")) %>%
  mutate(source = id, target = "Hypertension", interaction = "leads_to") %>%
  select(source, target, interaction)

# Orphan pollutants -> disease
edges_orphan <- nodes_cy %>%
  filter(Type %in% names(user_colors)) %>%
  filter(!id %in% linked_pollutants) %>%
  mutate(source = id, target = "Hypertension", interaction = "direct_to_disease") %>%
  select(source, target, interaction)

edges_cy <- bind_rows(edges_tmp, edges_bio, edges_orphan) %>% distinct()

# Keep only edges whose endpoints are in the node table
edges_cy <- edges_cy %>%
  filter(source %in% nodes_cy$id & target %in% nodes_cy$id)

cat("\nEdges after filtering:", nrow(edges_cy), "\n")
cat("Protein edges:", nrow(edges_gene), "\n")

# -----------------------------------------------------------------------------
# Data cleaning and network creation
# -----------------------------------------------------------------------------
nodes_cy <- nodes_cy %>%
  mutate(across(c(id, Type, Label, Is_Final), as.character),
         Base_Imp = as.numeric(Base_Imp),
         Plot_Size = as.numeric(Plot_Size))

edges_cy <- edges_cy %>%
  mutate(across(c(source, target, interaction), as.character)) %>%
  filter(!is.na(source), !is.na(target))

net_name <- "Hypertension Core Network (SG only)"
collection <- "AOP Framework"
existing <- tryCatch(getNetworkList(), error = function(e) list())
if (is.list(existing) && net_name %in% names(existing)) {
  cat("Deleting existing network with the same name...\n")
  deleteNetwork(existing[net_name])
  Sys.sleep(0.5)
}

tryCatch({
  createNetworkFromDataFrames(
    nodes = nodes_cy,
    edges = edges_cy,
    title = net_name,
    collection = collection
  )
}, error = function(e) {
  stop("Network creation failed. Error: ", e$message)
})

# -----------------------------------------------------------------------------
# Visual style
# -----------------------------------------------------------------------------
style_name <- "Hypertension_Core_Style"
try(deleteVisualStyle(style_name), silent = TRUE)

defaults <- list(
  NODE_SHAPE = "ELLIPSE",
  NODE_LABEL_FONT_SIZE = 14,
  NODE_LABEL_COLOR = "#000000",
  EDGE_TRANSPARENCY = 255,
  EDGE_WIDTH = 2.0,
  EDGE_TARGET_ARROW_SHAPE = "NONE",
  NODE_BORDER_PAINT = "#000000",
  NODE_BORDER_WIDTH = 8.0,
  EDGE_PAINT = "#BEBEBE"
)
createVisualStyle(style_name, defaults, mapping = list())
setVisualStyle(style_name)

all_types <- c(names(user_colors), "GO", "Pathway", "Protein", "Disease", "Pollutant")
all_colors <- c(user_colors,
                "GO"        = "#A9A7D6",
                "Pathway"   = "#66CDAA",
                "Protein"   = "#F8D582",
                "Disease"   = "#E64B35",
                "Pollutant" = "#808080")
setNodeColorMapping(
  table.column = 'Type',
  table.column.values = all_types,
  colors = all_colors,
  mapping.type = "d",
  style.name = style_name
)

min_size <- min(nodes_cy$Plot_Size[nodes_cy$Type != "Disease"], na.rm = TRUE)
max_size <- max(nodes_cy$Plot_Size[nodes_cy$Type != "Disease"], na.rm = TRUE)
setNodeSizeMapping(
  table.column = 'Plot_Size',
  table.column.values = c(min_size, max_size),
  sizes = c(30, 120),
  mapping.type = "c",
  style.name = style_name
)

setNodeLabelMapping('Label', style.name = style_name)

layoutNetwork('force-directed')

# -----------------------------------------------------------------------------
# Legend
# -----------------------------------------------------------------------------
legend_df <- data.frame(Type = factor(all_types, levels = all_types), x = 1, y = 1)
p_legend <- ggplot(legend_df, aes(x = x, y = y, fill = Type)) +
  geom_point(shape = 21, size = 6, color = "#555555", stroke = 1) +
  scale_fill_manual(values = all_colors, name = "Node Fill Color") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold", margin = margin(b = 10)),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm")
  ) +
  guides(fill = guide_legend(ncol = 1))
leg <- cowplot::get_legend(p_legend)
ggsave("Core_Network_Legend.pdf", plot = leg, width = 4, height = 8, bg = "white")

cat("\nCore network built.\n")
cat("Nodes (SG core features only):", nrow(nodes_cy),
    ", edges:", nrow(edges_cy), "\n")
