# =============================================================================
# ROC curves with DeLong tests (five SG models) + feature importance donuts
# =============================================================================

suppressMessages({
  library(mlr3verse)
  library(data.table)
  library(pROC)
  library(ggplot2)
  library(ggsci)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# Input (load before running): saved model results
#   Eno-D  <- readRDS("BD/eSet_BD.rds")$model
#   E-D    <- readRDS("ED/eSet_ED.rds")$ED_model
#   E-GO-D <- readRDS("GO/eSet_result.rds")$sg
#   E-P-D  <- readRDS("gene153/eSet_result.rds")$sg
#   E-Pa-D <- readRDS("pathway_v4/eSet_result.rds")$sg
# -----------------------------------------------------------------------------

get_pred <- function(x, nm) {
  rr <- x[[nm]]$bmr_final$resample_result(1)
  dt <- as.data.table(rr$prediction())
  list(truth = dt$truth, prob = dt$prob.1)
}

models <- list(
  `Eno-D`  = get_pred(EnoD, "model"),
  `E-D`    = get_pred(ED, "ED_model"),
  `E-GO-D` = get_pred(EGOD, "sg"),
  `E-P-D`  = get_pred(EPD, "sg"),
  `E-Pa-D` = get_pred(EPaD, "sg")
)

rocs <- lapply(models, function(p) roc(p$truth, p$prob, quiet = TRUE))

labels <- sapply(names(rocs), function(nm) {
  c <- ci.auc(rocs[[nm]], method = "delong")
  sprintf("%s (AUC = %.3f, 95%% CI: %.3f-%.3f)", nm, as.numeric(rocs[[nm]]$auc), c[1], c[3])
})

roc_df <- do.call(rbind, lapply(names(rocs), function(nm) {
  data.frame(FPR = 1 - rocs[[nm]]$specificities, TPR = rocs[[nm]]$sensitivities,
             Model = labels[nm], stringsAsFactors = FALSE)
}))
roc_df$Model <- factor(roc_df$Model, levels = labels)

# Pairwise DeLong tests vs E-D
p_delong <- sapply(c("E-GO-D", "E-P-D", "E-Pa-D"), function(nm) {
  t <- roc.test(rocs[["E-D"]], rocs[[nm]], method = "delong", paired = TRUE)
  t$p.value
})

subtitle <- sprintf("DeLong test vs E-D: GO P=%.3f | P(153) P=%.3f | Pa P=%.3f",
                    p_delong["E-GO-D"], p_delong["E-P-D"], p_delong["E-Pa-D"])

# -----------------------------------------------------------------------------
# ROC figure
# -----------------------------------------------------------------------------
p1 <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60", linewidth = 0.8) +
  geom_path(linewidth = 1.1, alpha = 0.92) +
  scale_color_jco() +
  labs(title = "Predictive Performance of the Five Models",
       subtitle = subtitle,
       x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)") +
  theme_bw() +
  theme(legend.position = c(0.68, 0.22),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        plot.subtitle = element_text(size = 9, color = "gray30"))

ggsave("Figure3a_ROC.png", p1, width = 8, height = 7, dpi = 600)
ggsave("Figure3a_ROC.pdf", p1, width = 8, height = 7)

# -----------------------------------------------------------------------------
# Feature importance donut charts
# -----------------------------------------------------------------------------
groups_fine <- list(
  Sulfonamides     = c("SFM", "4AN6C2B", "SFSX", "TMTH", "SFTHI"),
  Fluoroquinolones = c("NOR", "ENOR", "CIP", "OFLX"),
  Tetracyclines    = c("TET", "OTC", "AM", "DOXY"),
  Neonicotinoids   = c("IMI", "NIT", "ACE", "THM", "CLO", "THI", "IPP", "CYC",
                       "DIN", "IMID", "DM-THM", "DN-IMI", "5-OH-IMI", "DM-ACE",
                       "IMI-olefin"),
  Organophosphates = c("DZN", "CPS-O", "3,5,6-TCP", "DDVP", "MALT"),
  Carbamates       = c("3OH-CBF"),
  Pyrethroids      = c("3-PBA", "4-F-3-PBA", "2-MB-3CA", "trans-DCCA", "cis-DCCA"),
  Herbicides       = c("ATZ", "ATZ-DE-DSI-2HD", "CBDZ", "DMM"),
  Phenylpyrazoles  = c("FP"),
  Bisphenols       = c("BPF", "BPS", "BPC", "BPE", "BPP"),
  OH_PAHs          = c("2-OHN", "1-OHN", "3-OHF", "2-OHF", "2-OHP", "1-OHP",
                       "4-OHP", "1-OHPYR"),
  Metals           = c("V", "Fe", "Co", "Cu", "Zn", "As", "Se", "Sr", "Mo",
                       "Cd", "Te", "Tl", "Pb")
)

user_colors <- c(
  "Sulfonamides"     = "#C0392B", "Fluoroquinolones" = "#E74C3C",
  "Tetracyclines"    = "#E67E22", "Neonicotinoids" = "#F39C12",
  "Organophosphates" = "#F1C40F", "Carbamates" = "#D4AC0D",
  "Pyrethroids"      = "#82E0AA", "Herbicides" = "#58D68D",
  "Phenylpyrazoles"  = "#2ECC71", "Bisphenols" = "#5DADE2",
  "OH_PAHs"          = "#3498DB", "Metals" = "#1A5276",
  "GO" = "#A9A7D6", "Protein" = "#F8D582", "Pathway" = "#66CDAA",
  "Baseline" = "#A9A9A9", "Clinical" = "#FF8C00",
  "Base model" = "#FF4500", "Other" = "#DDDDDD"
)

pollutant_map <- data.frame(
  RawName = unlist(groups_fine),
  Group = rep(names(groups_fine), lengths(groups_fine)),
  stringsAsFactors = FALSE
)

classify <- function(f) {
  case_when(
    grepl("^pred\\.1\\.", f) ~ "Base model",
    grepl("^Baseline_", f) ~ "Baseline",
    grepl("^Clinical_", f) ~ "Clinical",
    grepl("^Exposome_GO_", f) ~ "GO",
    grepl("^Exposome_Gene_", f) ~ "Protein",
    grepl("^Exposome_Pathway_", f) ~ "Pathway",
    grepl("^Exposome_", f) ~ {
      raw <- sub("^Exposome_", "", f)
      g <- pollutant_map$Group[match(raw, pollutant_map$RawName)]
      ifelse(is.na(g), "Other", g)
    },
    TRUE ~ "Other"
  )
}

# Input: imp lists per model
# imp_list <- list(`E-GO-D` = EGOD$sg$imp_df, `E-P-D` = EPD$sg$imp_df,
#                  `E-Pa-D` = EPaD$sg$imp_df)

for (nm in names(imp_list)) {
  imp <- imp_list[[nm]] %>% filter(!is.na(Importance)) %>%
    mutate(Type = classify(Features),
           Label = case_when(
             Type %in% c("GO", "Protein", "Pathway") ~ sub("^Exposome_(GO|Gene|Pathway)_", "", Features),
             Type == "Base model" ~ sub("^pred\\.1\\.", "", Features),
             TRUE ~ sub("^Exposome_|^Baseline_|^Clinical_", "", Features)
           ))
  imp <- imp %>% arrange(desc(Importance_pct))
  imp <- imp %>% mutate(Type = ifelse(Importance_pct < 0.02, "Other", Type))
  imp$Type <- factor(imp$Type, levels = names(user_colors))

  p <- ggplot(imp, aes(x = 2, y = Importance_pct, fill = Type)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.3) +
    coord_polar(theta = "y", start = 0) +
    xlim(0.5, 2.5) +
    scale_fill_manual(values = user_colors, drop = FALSE, name = "Feature type") +
    labs(title = sprintf("Feature Importance - %s", nm)) +
    theme_void() +
    theme(legend.position = "right",
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

  ggsave(sprintf("Figure3b_Donut_%s.png", nm), p, width = 8, height = 6, dpi = 600)
  ggsave(sprintf("Figure3b_Donut_%s.pdf", nm), p, width = 8, height = 6)
}

cat("AUCs:\n")
print(labels)
cat(subtitle, "\ndone\n")
