# =============================================================================
# Weighted quantile sum (WQS) regression: positive and negative directions
# Outcome: hypertension (binary)
# Top-20 weight bar plots by chemical class (1:2 aspect ratio)
# =============================================================================

rm(list = ls())

if (!require("gWQS")) install.packages("gWQS")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("patchwork")) install.packages("patchwork")
if (!require("dplyr")) install.packages("dplyr")

library(gWQS)
library(ggplot2)
library(patchwork)
library(dplyr)

# -----------------------------------------------------------------------------
# Input (load before running)
# -----------------------------------------------------------------------------
# df <- read.xlsx(...)

# -----------------------------------------------------------------------------
# Variables
# -----------------------------------------------------------------------------
pollutants <- c(
  "SFM", "4AN6C2B", "SFSX", "TMTH", "SFTHI",
  "NOR", "ENOR", "CIP", "OFLX", "TET",
  "OTC", "AM", "DOXY", "IMI", "NIT",
  "ACE", "THM", "CLO", "THI", "IPP",
  "CYC", "DIN", "IMID", "DM-THM", "DN-IMI",
  "5-OH-IMI", "DM-ACE", "FP", "DMM", "CBDZ",
  "CPS-O", "DZN", "DDVP", "MALT", "3OH-CBF",
  "ATZ", "ATZ-DE-DSI-2HD", "BPF", "BPS",
  "BPC", "BPE", "BPP", "2-OHN", "1-OHN",
  "3-OHF", "2-OHF", "2-OHP", "1-OHP", "4-OHP",
  "1-OHPYR", "3-PBA", "4-F-3-PBA", "2-MB-3CA",
  "trans-DCCA", "cis-DCCA", "IMI-olefin", "3,5,6-TCP",
  "V", "Fe", "Co", "Cu", "Zn", "As", "Se",
  "Sr", "Mo", "Cd", "Te", "Tl", "Pb"
)

covariates_cat <- c("Gender", "Exercisefrequency", "Eatinghabits", "Smoke", "Drink")
covariates_cont <- c("Leftfastingbloodglucose2mmol", "BMI", "Age2")
all_covariates <- c(covariates_cat, covariates_cont)
outcome <- "Hypertension"

pollutant_groups <- list(
  Sulfonamides     = c("SFM", "4AN6C2B", "SFSX", "TMTH", "SFTHI"),
  Fluoroquinolones = c("NOR", "ENOR", "CIP", "OFLX"),
  Tetracyclines    = c("TET", "OTC", "AM", "DOXY"),
  Neonicotinoids   = c("IMI", "NIT", "ACE", "THM", "CLO", "THI", "IPP", "CYC", "DIN", "IMID",
                       "DM-THM", "DN-IMI", "5-OH-IMI", "DM-ACE", "IMI-olefin"),
  Organophosphates = c("DZN", "CPS-O", "3,5,6-TCP", "DDVP", "MALT"),
  Carbamates       = c("3OH-CBF"),
  Pyrethroids      = c("3-PBA", "4-F-3-PBA", "2-MB-3CA", "trans-DCCA", "cis-DCCA"),
  Herbicides       = c("ATZ", "ATZ-DE-DSI-2HD", "CBDZ", "DMM"),
  Phenylpyrazoles  = c("FP"),
  Bisphenols       = c("BPF", "BPS", "BPC", "BPE", "BPP"),
  OHPAHs           = c("2-OHN", "1-OHN", "3-OHF", "2-OHF", "2-OHP", "1-OHP", "4-OHP", "1-OHPYR"),
  Metals           = c("V", "Fe", "Co", "Cu", "Zn", "As", "Se", "Sr", "Mo", "Cd", "Te", "Tl", "Pb")
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
  "OHPAHs"           = "#3498DB",
  "Metals"           = "#1A5276"
)

get_group <- function(name) {
  for (g in names(pollutant_groups)) {
    if (name %in% pollutant_groups[[g]]) return(g)
  }
  return(NA)
}

# -----------------------------------------------------------------------------
# Merge highly collinear variables
# -----------------------------------------------------------------------------
high_corr_pairs <- list(
  DCCA = c("trans-DCCA", "cis-DCCA"),
  TETDOXY = c("TET", "DOXY")
)

for (pair in names(high_corr_pairs)) {
  vars <- high_corr_pairs[[pair]]
  if (all(vars %in% names(df))) {
    df[[pair]] <- rowMeans(df[, vars], na.rm = TRUE)
    cat(paste("Merged variable:", pair, "= mean of", paste(vars, collapse = " + "), "\n"))
    pollutants <- c(pollutants, pair)
    pollutants <- pollutants[!pollutants %in% vars]
  } else {
    warning(paste("Original variables for", pair, "missing; merge skipped."))
  }
}

pollutants <- unique(pollutants)

pollutants_exist <- pollutants[pollutants %in% names(df)]
if (length(pollutants_exist) == 0) stop("No pollutant variables found; check names.")
cat("Pollutants used (including merged):", length(pollutants_exist), "\n")

if ("DCCA" %in% pollutants_exist) {
  pollutant_groups[["Pyrethroids"]] <- c("3-PBA", "4-F-3-PBA", "2-MB-3CA", "DCCA")
}
if ("TETDOXY" %in% pollutants_exist) {
  pollutant_groups[["Tetracyclines"]] <- c("OTC", "AM", "TETDOXY")
}

all_vars <- c(outcome, pollutants_exist, all_covariates)
missing_vars <- all_vars[!all_vars %in% names(df)]
if (length(missing_vars) > 0) {
  warning("Variables not found: ", paste(missing_vars, collapse = ", "))
  pollutants_exist <- setdiff(pollutants_exist, missing_vars)
  all_covariates <- setdiff(all_covariates, missing_vars)
}

df_model <- na.omit(df[, c(outcome, pollutants_exist, all_covariates)])
cat("Sample size after removing missing:", nrow(df_model), "\n")

for (v in covariates_cat) {
  if (v %in% names(df_model) && !is.factor(df_model[[v]])) {
    df_model[[v]] <- as.factor(df_model[[v]])
  }
}

if (!is.numeric(df_model$Hypertension)) {
  df_model$Hypertension <- as.numeric(df_model$Hypertension)
}
unique_vals <- unique(df_model$Hypertension)
if (!all(unique_vals %in% c(0, 1))) {
  stop("Outcome Hypertension is not binary 0/1; check the original coding.")
}

# -----------------------------------------------------------------------------
# WQS positive direction (4 quantiles)
# -----------------------------------------------------------------------------
formula_wqs <- as.formula(paste(outcome, "~ wqs +", paste(all_covariates, collapse = " + ")))

set.seed(123)
system.time({
  wqs_pos <- gwqs(formula_wqs,
                  mix_name = pollutants_exist,
                  data = df_model,
                  q = 4,
                  validation = 0.4,
                  b = 500,
                  b1_pos = TRUE,
                  b_constr = TRUE,
                  family = "binomial",
                  seed = 123,
                  plots = FALSE,
                  tables = FALSE)
})

weights_pos <- wqs_pos$final_weights
weights_pos$mix_name_clean <- weights_pos$mix_name
weights_pos$group <- sapply(weights_pos$mix_name_clean, get_group)
if (any(is.na(weights_pos$group))) {
  weights_pos$group[is.na(weights_pos$group)] <- "Other"
  user_colors <- c(user_colors, Other = "#999999")
}

wqs_pos_summary <- gwqs_summary_tab(wqs_pos)

# -----------------------------------------------------------------------------
# WQS negative direction (4 quantiles)
# -----------------------------------------------------------------------------
set.seed(456)
system.time({
  wqs_neg <- gwqs(formula_wqs,
                  mix_name = pollutants_exist,
                  data = df_model,
                  q = 4,
                  validation = 0.4,
                  b = 500,
                  b1_pos = FALSE,
                  b_constr = TRUE,
                  family = "binomial",
                  seed = 456,
                  plots = FALSE,
                  tables = FALSE)
})

weights_neg <- wqs_neg$final_weights
weights_neg$mix_name_clean <- weights_neg$mix_name
weights_neg$group <- sapply(weights_neg$mix_name_clean, get_group)
if (any(is.na(weights_neg$group))) {
  weights_neg$group[is.na(weights_neg$group)] <- "Other"
}

wqs_neg_summary <- gwqs_summary_tab(wqs_neg)

# -----------------------------------------------------------------------------
# Top-20 weights
# -----------------------------------------------------------------------------
top_n_weights <- function(weights_df, n = 20) {
  weights_df %>%
    arrange(desc(mean_weight)) %>%
    slice_head(n = n) %>%
    arrange(mean_weight)
}

weights_pos_top <- top_n_weights(weights_pos, n = 20)
weights_neg_top <- top_n_weights(weights_neg, n = 20)

cat("Positive direction variables plotted:", nrow(weights_pos_top), "\n")
cat("Negative direction variables plotted:", nrow(weights_neg_top), "\n")

# -----------------------------------------------------------------------------
# Bar plots
# -----------------------------------------------------------------------------
theme_bar <- function() {
  theme_bw() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
}

p1 <- ggplot(weights_pos_top, aes(x = reorder(mix_name_clean, mean_weight),
                                  y = mean_weight, fill = group)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = user_colors, name = "Chemical Class") +
  coord_flip() +
  labs(title = "WQS Positive Direction (Top 20 Weights)", x = "Weight", y = NULL) +
  theme_bar()

p2 <- ggplot(weights_neg_top, aes(x = reorder(mix_name_clean, mean_weight),
                                  y = mean_weight, fill = group)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = user_colors, name = "Chemical Class") +
  coord_flip() +
  labs(title = "WQS Negative Direction (Top 20 Weights)", x = "Weight", y = NULL) +
  theme_bar()

combined_plot <- p1 + p2 + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave("WQS_Combined_Weights_Top20.pdf", combined_plot,
       width = 12, height = 16, dpi = 300, limitsize = FALSE)

ggsave("WQS_positive_weights_top20.pdf", p1,
       width = 6, height = 16, dpi = 300, limitsize = FALSE)
ggsave("WQS_negative_weights_top20.pdf", p2,
       width = 6, height = 16, dpi = 300, limitsize = FALSE)

# -----------------------------------------------------------------------------
# WQS index summaries
# -----------------------------------------------------------------------------
cat("\n===== WQS positive model summary =====\n")
print(wqs_pos_summary)

cat("\n===== WQS negative model summary =====\n")
print(wqs_neg_summary)

or_pos <- exp(coef(wqs_pos)[["wqs"]])
ci_pos <- exp(confint(wqs_pos, parm = "wqs"))
cat("\n===== Positive WQS index =====\n")
cat("OR (per quantile increase):", round(or_pos, 3), "\n")
cat("95% CI:", round(ci_pos[1], 3), "-", round(ci_pos[2], 3), "\n")
cat("P-value:", round(summary(wqs_pos)$coefficients["wqs", "Pr(>|z|)"], 4), "\n")

or_neg <- exp(coef(wqs_neg)[["wqs"]])
ci_neg <- exp(confint(wqs_neg, parm = "wqs"))
cat("\n===== Negative WQS index =====\n")
cat("OR (per quantile increase):", round(or_neg, 3), "\n")
cat("95% CI:", round(ci_neg[1], 3), "-", round(ci_neg[2], 3), "\n")
cat("P-value:", round(summary(wqs_neg)$coefficients["wqs", "Pr(>|z|)"], 4), "\n")
