# =============================================================================
# Biological feature matrix construction (ExposomeX platform)
# Pollutant x Gene (EGeD), pollutant x GO (EGoD) and pollutant x pathway
# (EPaD) matrices; individual perturbation matrices = exposure %*% matrix
# =============================================================================

rm(list = ls())

pacman::p_load("tidyverse", "readxl", "openxlsx")

# -----------------------------------------------------------------------------
# Input (load before running):
#   tidy_df: ExposomeX mapping table (columns: abbreviation, exid)
#   eged_df / egod_df / epad_df: edge tables exported from ExposomeX
#   (columns: Source, Target, "Source group", "Target group", ...)
#   exposure_df: exposure data (pollutant columns with creatinine-adjustment
#   suffix "_jgjz")
# -----------------------------------------------------------------------------
# tidy_df <- read_excel("tidy_data.xlsx")
# eged_df <- read_csv("EGeD_edges.csv")
# egod_df <- read_csv("EGoD_edges.csv")
# epad_df <- read_csv("EPaD_edges.csv")
# exposure_df <- read_excel("hypertension.xlsx")

# -----------------------------------------------------------------------------
# Pollutant abbreviation -> ExposomeX exposure ID
# -----------------------------------------------------------------------------
short_to_ex <- tidy_df %>%
  dplyr::select(abbreviation, exid) %>%
  dplyr::filter(!is.na(abbreviation), !is.na(exid),
                abbreviation != "nan", exid != "nan") %>%
  deframe()

# Manual metal mapping
metal_mapping <- c(
  "V"  = "EX:E000022460",
  "Fe" = "EX:E000022401",
  "Co" = "EX:E000096046",
  "Cu" = "EX:E000022448",
  "Zn" = "EX:E000022464",
  "As" = "EX:E004815849",
  "Se" = "EX:E005270968",
  "Sr" = "EX:E004815779",
  "Mo" = "EX:E000022407",
  "Cd" = "EX:E000022444",
  "Te" = "EX:E005271068",
  "Tl" = "EX:E004815817",
  "Pb" = "EX:E004811303"
)
short_to_ex <- c(short_to_ex, metal_mapping)
print(paste("Mapping entries (with metals):", length(short_to_ex)))

# -----------------------------------------------------------------------------
# Exposure matrix
# -----------------------------------------------------------------------------
exposure_cols <- grep("_jgjz", names(exposure_df), value = TRUE)
print(paste("Exposure columns:", length(exposure_cols)))

valid_cols <- c()
valid_ex <- c()
valid_shorts <- c()

for (col in exposure_cols) {
  short <- stringr::str_replace(col, "_\\d*_jgjz$", "")
  short <- stringr::str_replace(short, "_jgjz$", "")

  if (short %in% names(short_to_ex)) {
    valid_cols <- c(valid_cols, col)
    valid_ex <- c(valid_ex, short_to_ex[short])
    valid_shorts <- c(valid_shorts, short)
  } else {
    cat("Unmapped column:", col, "-> abbreviation:", short, "\n")
  }
}

exposure_mat <- exposure_df %>%
  dplyr::select(all_of(valid_cols)) %>%
  as.matrix()
colnames(exposure_mat) <- valid_ex

print(paste("Exposure matrix:", nrow(exposure_mat), "rows,",
            ncol(exposure_mat), "columns"))
print(paste("Valid pollutants:", length(valid_ex)))

metal_cols <- c("V", "Fe", "Co", "Cu", "Zn", "As", "Se", "Sr", "Mo", "Cd",
                "Te", "Tl", "Pb")
present_metals <- metal_cols[metal_cols %in% valid_shorts]
print(paste("Metals mapped:", paste(present_metals, collapse = ", ")))

# -----------------------------------------------------------------------------
# Pollutant x Gene matrix (EGeD)
# -----------------------------------------------------------------------------
gene_matrix <- eged_df %>%
  dplyr::select(Source, Target) %>%
  dplyr::filter(Source %in% valid_ex) %>%
  dplyr::mutate(value = 1) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(
    names_from = Target,
    names_sep = "_",
    values_from = value,
    values_fill = 0
  ) %>%
  dplyr::rename(name = Source)

print(paste("Pollutant x gene matrix:", nrow(gene_matrix), "rows,",
            ncol(gene_matrix), "columns"))

# -----------------------------------------------------------------------------
# Pollutant x GO matrix (EGoD)
# -----------------------------------------------------------------------------
go_matrix <- egod_df %>%
  dplyr::select(Source, Target) %>%
  dplyr::filter(Source %in% valid_ex) %>%
  dplyr::mutate(value = 1) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(
    names_from = Target,
    names_sep = "_",
    values_from = value,
    values_fill = 0
  ) %>%
  dplyr::rename(name = Source)

print(paste("Pollutant x GO matrix:", nrow(go_matrix), "rows,",
            ncol(go_matrix), "columns"))

# -----------------------------------------------------------------------------
# Pollutant x pathway matrix (EPaD)
# -----------------------------------------------------------------------------
pathway_matrix <- epad_df %>%
  dplyr::select(Source, Target) %>%
  dplyr::filter(Source %in% valid_ex) %>%
  dplyr::mutate(value = 1) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(
    names_from = Target,
    names_sep = "_",
    values_from = value,
    values_fill = 0
  ) %>%
  dplyr::rename(name = Source)

print(paste("Pollutant x pathway matrix:", nrow(pathway_matrix), "rows,",
            ncol(pathway_matrix), "columns"))

# -----------------------------------------------------------------------------
# Individual perturbation matrices (exposure %*% biological matrix)
# -----------------------------------------------------------------------------
if (exists("gene_matrix") && nrow(gene_matrix) > 0) {
  gene_mat <- gene_matrix %>%
    dplyr::select(-name) %>%
    as.matrix()
  rownames(gene_mat) <- gene_matrix$name

  common_pollutants_gene <- intersect(colnames(exposure_mat), rownames(gene_mat))
  print(paste("Common pollutants (gene):", length(common_pollutants_gene)))

  if (length(common_pollutants_gene) > 0) {
    exposure_sub <- as.matrix(exposure_mat[, common_pollutants_gene, drop = FALSE])
    gene_sub <- as.matrix(gene_mat[common_pollutants_gene, , drop = FALSE])
    individual_gene <- exposure_sub %*% gene_sub
    colnames(individual_gene) <- colnames(gene_sub)
  }
}

if (exists("go_matrix") && nrow(go_matrix) > 0) {
  go_mat <- go_matrix %>%
    dplyr::select(-name) %>%
    as.matrix()
  rownames(go_mat) <- go_matrix$name

  common_pollutants_go <- intersect(colnames(exposure_mat), rownames(go_mat))
  print(paste("Common pollutants (GO):", length(common_pollutants_go)))

  if (length(common_pollutants_go) > 0) {
    exposure_sub <- as.matrix(exposure_mat[, common_pollutants_go, drop = FALSE])
    go_sub <- as.matrix(go_mat[common_pollutants_go, , drop = FALSE])
    individual_go <- exposure_sub %*% go_sub
    colnames(individual_go) <- colnames(go_sub)
  }
}

if (exists("pathway_matrix") && nrow(pathway_matrix) > 0) {
  pathway_mat <- pathway_matrix %>%
    dplyr::select(-name) %>%
    as.matrix()
  rownames(pathway_mat) <- pathway_matrix$name

  common_pollutants_path <- intersect(colnames(exposure_mat), rownames(pathway_mat))
  print(paste("Common pollutants (pathway):", length(common_pollutants_path)))

  if (length(common_pollutants_path) > 0) {
    exposure_sub <- as.matrix(exposure_mat[, common_pollutants_path, drop = FALSE])
    pathway_sub <- as.matrix(pathway_mat[common_pollutants_path, , drop = FALSE])
    individual_pathway <- exposure_sub %*% pathway_sub
    colnames(individual_pathway) <- colnames(pathway_sub)
  }
}

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
results <- list(
  hypertension = exposure_df$Hypertension,
  individual_gene = if (exists("individual_gene")) individual_gene else NULL,
  individual_go = if (exists("individual_go")) individual_go else NULL,
  individual_pathway = if (exists("individual_pathway")) individual_pathway else NULL,
  pollutant_gene_matrix = if (exists("gene_matrix")) gene_matrix else NULL,
  pollutant_go_matrix = if (exists("go_matrix")) go_matrix else NULL,
  pollutant_pathway_matrix = if (exists("pathway_matrix")) pathway_matrix else NULL
)

cat("\nPollutants:", ncol(exposure_mat), "\n")
if (exists("individual_gene") && !is.null(individual_gene)) {
  cat("Genes:", ncol(individual_gene), "\n")
} else {
  cat("Genes: 0\n")
}
if (exists("individual_go") && !is.null(individual_go)) {
  cat("GO terms:", ncol(individual_go), "\n")
} else {
  cat("GO terms: 0\n")
}
if (exists("individual_pathway") && !is.null(individual_pathway)) {
  cat("Pathways:", ncol(individual_pathway), "\n")
} else {
  cat("Pathways: 0\n")
}
cat("Samples:", length(exposure_df$Hypertension), "\n")
cat("Hypertension cases:", sum(exposure_df$Hypertension), "\n")
cat("Controls:", sum(1 - exposure_df$Hypertension), "\n")
