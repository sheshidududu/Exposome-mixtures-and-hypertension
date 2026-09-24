# =============================================================================
# Spearman correlation network of 70 pollutants (12 classes) and hypertension,
# exported to Cytoscape; intra-group / all-pair correlation count bar charts
# =============================================================================

library(tidyverse)
library(readxl)
library(writexl)
library(rstatix)
library(purrr)
library(RCy3)
library(ggplot2)
library(tidyr)
library(tidytext)   # reorder_within / scale_y_reordered

# -----------------------------------------------------------------------------
# Input (load before running)
# -----------------------------------------------------------------------------
# data <- read_excel("hypertension_imputed.xlsx")
if (!"Hypertension" %in% names(data)) {
  stop("Hypertension variable missing; check column names")
}
data$Hypertension <- as.numeric(data$Hypertension)

# -----------------------------------------------------------------------------
# Exposure variable groups (12 classes)
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

varlist <- unlist(groups_fine) %>% unique()
missing_vars <- varlist[!varlist %in% names(data)]
if (length(missing_vars) > 0) {
  warning("Variables not found: ", paste(missing_vars, collapse = ", "))
  varlist <- varlist[varlist %in% names(data)]
}

var_to_group <- tibble(
  variable = varlist,
  group = map_chr(varlist, function(x) {
    grp <- names(groups_fine)[map_lgl(groups_fine, ~ x %in% .x)]
    if (length(grp) == 0) NA_character_ else grp[1]
  })
)

# -----------------------------------------------------------------------------
# Spearman correlations between exposures
# -----------------------------------------------------------------------------
comb <- combn(varlist, 2) %>% t() %>% as_tibble() %>% setNames(c("var1", "var2"))

cor_func <- function(v1, v2) {
  cor_test(data, vars = v1, vars2 = v2, method = "spearman")
}
cor_result <- map2_dfr(comb$var1, comb$var2, cor_func) %>%
  mutate(
    group1 = var_to_group$group[match(var1, var_to_group$variable)],
    group2 = var_to_group$group[match(var2, var_to_group$variable)]
  )

# -----------------------------------------------------------------------------
# Exposure-outcome correlation (node size)
# -----------------------------------------------------------------------------
node_size <- map_dfr(varlist, function(x) {
  cor_test(data, vars = "Hypertension", vars2 = x, method = "spearman")
}) %>%
  mutate(MinusLogP = -log10(p + 1e-8), variable = var2) %>%
  select(variable, cor, p, MinusLogP)

# -----------------------------------------------------------------------------
# Node and edge tables
# -----------------------------------------------------------------------------
nodes <- var_to_group %>%
  left_join(node_size, by = "variable") %>%
  mutate(id = variable, label = variable) %>%
  select(id, label, group, cor, p, MinusLogP)

edges <- cor_result %>%
  filter(p < 0.05) %>%
  mutate(source = var1, target = var2, interaction = "association") %>%
  select(source, target, interaction, cor, p) %>%
  distinct(source, target, .keep_all = TRUE) %>%
  left_join(cor_result %>% select(var1, var2, group1, group2),
            by = c("source" = "var1", "target" = "var2"))

cat("Nodes:", nrow(nodes), " Edges:", nrow(edges), "\n")

# -----------------------------------------------------------------------------
# Cytoscape
# -----------------------------------------------------------------------------
RCy3::cytoscapePing()
RCy3::closeSession(F)
RCy3::createNetworkFromDataFrames(nodes, edges,
                                  title = "Exposure_Network_Hypertension",
                                  collection = "EnvironmentalExposures")

style_name <- "Hypertension_Style"
defaults <- list(NODE_SHAPE = "ELLIPSE", NODE_SIZE = 20, EDGE_TRANSPARENCY = 225,
                 NODE_LABEL_POSITION = "S,W,c,0.00,0.00")

nodeLabels <- RCy3::mapVisualProperty('node label', 'label', 'p')
groups_unique <- unique(nodes$group)
nodeShape <- RCy3::mapVisualProperty('Node Shape', 'group', "d",
                                     groups_unique, rep("ELLIPSE", length(groups_unique)))

size_range <- range(nodes$MinusLogP, na.rm = TRUE)
if (diff(size_range) == 0) size_range <- c(0, 1)
nodeSize <- RCy3::mapVisualProperty('Node Size', 'MinusLogP', "c",
                                    size_range, c(30, 150))

nodeLabelFontSize <- RCy3::mapVisualProperty('Node Label Font Size', 'group', "d",
                                             groups_unique, rep(12, length(groups_unique)))
nodeWidth <- RCy3::mapVisualProperty('Node Width', 'group', "d",
                                     groups_unique, rep(30, length(groups_unique)))
nodeZLocation <- RCy3::mapVisualProperty('Node Z Location', "group", 'd',
                                         groups_unique, rep(0, length(groups_unique)))

RCy3::createVisualStyle(style_name, defaults,
                        list(nodeLabels, nodeShape, nodeSize,
                             nodeLabelFontSize, nodeWidth, nodeZLocation))

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

matched_colors <- user_colors[groups_unique]
matched_colors[is.na(matched_colors)] <- "#CCCCCC"

RCy3::setNodeColorMapping(mapping.type = 'd', table.column = 'group',
                          table.column.values = groups_unique,
                          colors = matched_colors, style.name = style_name)
RCy3::setVisualStyle(style_name)
RCy3::layoutNetwork("force-directed")
RCy3::fitContent()

# -----------------------------------------------------------------------------
# Correlation count bar charts
# -----------------------------------------------------------------------------
cor_intra_counts <- cor_result %>%
  filter(group1 == group2) %>%
  mutate(abs_cor = abs(cor),
         cor_bin = cut(abs_cor, breaks = c(-Inf, 0.3, 0.6, Inf),
                       labels = c("<0.3", "0.3-0.6", ">0.6"), include.lowest = TRUE)) %>%
  count(group = group1, cor_bin, name = "frequency") %>%
  complete(group, cor_bin, fill = list(frequency = 0))

cor_all_counts <- cor_result %>%
  mutate(abs_cor = abs(cor),
         cor_bin = cut(abs_cor, breaks = c(-Inf, 0.3, 0.6, Inf),
                       labels = c("<0.3", "0.3-0.6", ">0.6"), include.lowest = TRUE)) %>%
  select(group1, group2, cor_bin) %>%
  pivot_longer(cols = c(group1, group2), names_to = "var_group", values_to = "group") %>%
  count(group, cor_bin, name = "frequency") %>%
  complete(group, cor_bin, fill = list(frequency = 0))

plot_facet <- function(data, title, filename) {
  p <- ggplot(data, aes(x = frequency,
                        y = reorder_within(group, frequency, cor_bin),
                        fill = group)) +
    geom_col(color = "black", linewidth = 0.25) +
    facet_wrap(~ cor_bin, ncol = 3, scales = "free_y") +
    scale_y_reordered() +
    scale_fill_manual(values = user_colors) +
    labs(x = "Frequency of Correlated Pairs", y = "Pollutant Group",
         fill = "Pollutant Group", title = title) +
    theme_bw() +
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          axis.title = element_text(size = 13),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          strip.text = element_text(size = 12, face = "bold"),
          legend.position = "bottom",
          panel.spacing = unit(0.8, "cm"))
  ggsave(filename, plot = p, width = 16, height = 9, device = "pdf", bg = "white")
}

plot_facet(cor_intra_counts,
           title = "Correlation Counts (Intra-group)",
           filename = "Bar_IntraGroup_Facet.pdf")
plot_facet(cor_all_counts,
           title = "Correlation Counts (All Pairs)",
           filename = "Bar_AllPairs_Facet.pdf")

p_dodge <- ggplot(cor_intra_counts, aes(x = frequency, y = cor_bin, fill = group)) +
  geom_col(position = position_dodge(width = 0.9), color = "black",
           linewidth = 0.25, width = 0.8) +
  scale_fill_manual(values = user_colors) +
  labs(title = "Correlation Counts (Intra-group)", x = "Frequency",
       y = "|Correlation|", fill = "Group") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 13), legend.position = "right")
ggsave("Bar_IntraGroup_Dodge.pdf", plot = p_dodge,
       width = 14, height = 7, device = "pdf", bg = "white")

# -----------------------------------------------------------------------------
# Strongest intra-group pairs and full correlation table
# -----------------------------------------------------------------------------
intra_cor <- cor_result %>%
  filter(group1 == group2) %>%
  mutate(abs_cor = abs(cor), group = group1) %>%
  select(group, var1, var2, cor, p, abs_cor)

strongest_intra <- intra_cor %>%
  group_by(group) %>%
  arrange(desc(abs_cor), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(rank = "strongest") %>%
  select(Group = group, Variable1 = var1, Variable2 = var2,
         Correlation = cor, P_value = p, Abs_Correlation = abs_cor)

print(strongest_intra)

overall_cor <- cor_result %>%
  select(Variable1 = var1, Variable2 = var2,
         Correlation = cor, P_value = p,
         Group1 = group1, Group2 = group2) %>%
  arrange(desc(abs(Correlation)))

expo_outcome <- node_size %>%
  left_join(var_to_group, by = "variable") %>%
  select(Exposure = variable, Group = group,
         Correlation = cor, P_value = p, MinusLogP) %>%
  arrange(desc(MinusLogP))

n_total <- nrow(expo_outcome)
n_sig <- sum(expo_outcome$P_value < 0.05, na.rm = TRUE)
cat(sprintf("Significant exposure-outcome correlations: %d / %d (%.1f%%)\n",
            n_sig, n_total, n_sig / n_total * 100))
