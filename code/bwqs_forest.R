# =============================================================================
# Bayesian weighted quantile sum (BWQS) analysis + single-pollutant GLM
# with forest plot; 12 chemical classes
# Multi-pollutant groups -> BWQS; single-pollutant groups -> GLM
# Highly collinear pairs merged: trans-DCCA + cis-DCCA -> DCCA;
# TET + DOXY -> TETDOXY
# =============================================================================

rm(list = ls())

# Load packages
library(BWQS)
library(dplyr)
library(ggplot2)
library(tidyr)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# -----------------------------------------------------------------------------
# Input (load before running)
# -----------------------------------------------------------------------------
# data <- read.xlsx(...)

# Variable type conversion
factor_vars <- c("Gender", "Exercisefrequency", "Eatinghabits", "Smoke", "Drink")
for (v in factor_vars) {
  if (v %in% names(data)) data[[v]] <- as.factor(data[[v]])
}
data$Age2 <- as.numeric(data$Age2)
data$BMI <- as.numeric(data$BMI)
data$Leftfastingbloodglucose2mmol <- as.numeric(data$Leftfastingbloodglucose2mmol)

# -----------------------------------------------------------------------------
# Outcome and covariates
# -----------------------------------------------------------------------------
outcome <- "Hypertension"
covariates <- c("Gender", "Age2", "BMI", "Exercisefrequency",
                "Eatinghabits", "Smoke", "Drink", "Leftfastingbloodglucose2mmol")

# -----------------------------------------------------------------------------
# Pollutant groups (12 classes)
# -----------------------------------------------------------------------------
pollutant_groups <- list(
  Sulfonamides     = c("SFM", "4AN6C2B", "SFSX", "TMTH", "SFTHI"),
  Fluoroquinolones = c("NOR", "ENOR", "CIP", "OFLX"),
  Tetracyclines    = c("TET", "OTC", "AM", "DOXY"),
  Neonicotinoids   = c("IMI", "NIT", "ACE", "THM", "CLO", "THI", "IPP",
                       "CYC", "DIN", "IMID", "DM-THM", "DN-IMI",
                       "5-OH-IMI", "DM-ACE", "IMI-olefin"),
  Organophosphates = c("DZN", "CPS-O", "3,5,6-TCP", "DDVP", "MALT"),
  Carbamates       = c("3OH-CBF"),
  Pyrethroids      = c("3-PBA", "4-F-3-PBA", "2-MB-3CA", "trans-DCCA", "cis-DCCA"),
  Herbicides       = c("ATZ", "ATZ-DE-DSI-2HD", "CBDZ", "DMM"),
  Phenylpyrazoles  = c("FP"),
  Bisphenols       = c("BPF", "BPS", "BPC", "BPE", "BPP"),
  OHPAHs           = c("2-OHN", "1-OHN", "3-OHF", "2-OHF", "2-OHP", "1-OHP", "4-OHP", "1-OHPYR"),
  Metals           = c("V", "Fe", "Co", "Cu", "Zn", "As", "Se", "Sr", "Mo", "Cd", "Te", "Tl", "Pb")
)

group_order <- names(pollutant_groups)

# -----------------------------------------------------------------------------
# Color scheme
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# Merge highly collinear variables
# -----------------------------------------------------------------------------
high_corr_pairs <- list(
  DCCA    = c("trans-DCCA", "cis-DCCA"),
  TETDOXY = c("TET", "DOXY")
)

merged_map <- list()

for (pair in names(high_corr_pairs)) {
  vars <- high_corr_pairs[[pair]]
  if (all(vars %in% names(data))) {
    data[[pair]] <- rowMeans(data[, vars], na.rm = TRUE)
    cat(sprintf("Merged variable: %s = mean(%s)\n", pair, paste(vars, collapse = " + ")))
    merged_map[[pair]] <- vars
  }
}

for (new_var in names(merged_map)) {
  old_vars <- merged_map[[new_var]]
  for (grp in names(pollutant_groups)) {
    if (any(old_vars %in% pollutant_groups[[grp]])) {
      pollutant_groups[[grp]] <- setdiff(pollutant_groups[[grp]], old_vars)
      if (!new_var %in% pollutant_groups[[grp]]) {
        pollutant_groups[[grp]] <- c(pollutant_groups[[grp]], new_var)
      }
      break
    }
  }
}

single_pollutant_groups <- names(pollutant_groups)[sapply(pollutant_groups, length) == 1]
multi_pollutant_groups <- names(pollutant_groups)[sapply(pollutant_groups, length) >= 2]

cat("\nSingle-pollutant groups (GLM):", paste(single_pollutant_groups, collapse = ", "), "\n")
cat("Multi-pollutant groups (BWQS):", paste(multi_pollutant_groups, collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# BWQS formula and result containers
# -----------------------------------------------------------------------------
formula_bwqs <- as.formula(paste(outcome, "~", paste(covariates, collapse = " + ")))

bwqs_results <- data.frame(
  Group  = character(),
  Model  = character(),
  OR     = numeric(),
  LCI    = numeric(),
  UCI    = numeric(),
  Rhat   = numeric(),
  Pvalue = numeric(),
  stringsAsFactors = FALSE
)
weight_results <- data.frame()

# -----------------------------------------------------------------------------
# Multi-pollutant groups -> BWQS
# -----------------------------------------------------------------------------
for (grp in multi_pollutant_groups) {
  cat("\n========================================\n")
  cat(sprintf("BWQS - %s\n", grp))
  cat("========================================\n")

  mix_vars <- pollutant_groups[[grp]]
  vars_needed <- c(outcome, covariates, mix_vars)
  df_model <- na.omit(data[, vars_needed])

  fit <- bwqs(
    formula   = formula_bwqs,
    mix_name  = mix_vars,
    data      = df_model,
    q         = 4,
    chains    = 3,
    iter      = 10000,
    thin      = 10,
    seed      = 123,
    family    = "binomial"
  )

  beta <- fit$summary_fit["beta1", "mean"]
  lci  <- fit$summary_fit["beta1", "2.5%"]
  uci  <- fit$summary_fit["beta1", "97.5%"]
  Rhat <- fit$summary_fit["beta1", "Rhat"]

  OR  <- exp(beta)
  LCI <- exp(lci)
  UCI <- exp(uci)
  pseudo_p <- ifelse(LCI > 1 | UCI < 1, 0.01, 0.5)

  bwqs_results <- rbind(bwqs_results, data.frame(
    Group  = grp, Model = "BWQS",
    OR = round(OR, 4), LCI = round(LCI, 4), UCI = round(UCI, 4),
    Rhat = round(Rhat, 4), Pvalue = pseudo_p,
    stringsAsFactors = FALSE
  ))

  weight_rows <- grep("^W_", rownames(fit$summary_fit))
  tmp_weight <- data.frame(
    Group    = grp,
    Chemical = sub("^W_", "", rownames(fit$summary_fit)[weight_rows]),
    Weight   = fit$summary_fit[weight_rows, "mean"],
    stringsAsFactors = FALSE
  )
  weight_results <- rbind(weight_results, tmp_weight)

  cat(sprintf("  OR = %.3f (%.3f - %.3f), Rhat = %.4f\n", OR, LCI, UCI, Rhat))
}

# -----------------------------------------------------------------------------
# Single-pollutant groups -> GLM
# -----------------------------------------------------------------------------
for (grp in single_pollutant_groups) {
  cat("\n========================================\n")
  cat(sprintf("GLM - %s\n", grp))
  cat("========================================\n")

  poll <- pollutant_groups[[grp]][1]
  if (!poll %in% names(data)) {
    warning(sprintf("Variable %s not in data; group %s skipped.", poll, grp))
    next
  }

  poll_safe <- paste0("`", poll, "`")
  form_glm <- as.formula(paste(outcome, "~", poll_safe, "+", paste(covariates, collapse = " + ")))
  fit_glm <- glm(form_glm, data = data, family = binomial)

  coef_sum <- coef(summary(fit_glm))
  if (nrow(coef_sum) < 2) next

  beta  <- coef_sum[2, 1]
  se    <- coef_sum[2, 2]
  p_val <- coef_sum[2, 4]

  OR  <- exp(beta)
  LCI <- exp(beta - 1.96 * se)
  UCI <- exp(beta + 1.96 * se)

  bwqs_results <- rbind(bwqs_results, data.frame(
    Group  = grp, Model = "GLM",
    OR = round(OR, 4), LCI = round(LCI, 4), UCI = round(UCI, 4),
    Rhat = NA, Pvalue = p_val,
    stringsAsFactors = FALSE
  ))

  weight_results <- rbind(weight_results, data.frame(
    Group = grp, Chemical = poll, Weight = 1.0,
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  %s: OR = %.3f (%.3f - %.3f), P = %.4f\n", poll, OR, LCI, UCI, p_val))
}

# -----------------------------------------------------------------------------
# Result tables
# -----------------------------------------------------------------------------
bwqs_results$Group <- factor(bwqs_results$Group, levels = group_order)
bwqs_results <- bwqs_results %>% arrange(Group)

group_sizes <- sapply(pollutant_groups, length)

S13 <- weight_results %>%
  mutate(Significance_Threshold = round(1 / group_sizes[Group], 6)) %>%
  arrange(Group, desc(Weight)) %>%
  transmute(
    Chemical_Group         = Group,
    Pollutant              = Chemical,
    Posterior_Mean_Weight  = round(Weight, 6),
    Significance_Threshold = Significance_Threshold
  )

S12 <- bwqs_results %>%
  transmute(
    Chemical_Group     = Group,
    Model              = Model,
    OR                 = round(OR, 3),
    `95%_CI_Lower`     = round(LCI, 3),
    `95%_CI_Upper`     = round(UCI, 3),
    Pseudo_or_Actual_P = ifelse(Model == "BWQS",
                                as.character(Pvalue),
                                ifelse(Pvalue < 0.001, "<0.001",
                                       as.character(round(Pvalue, 4))))
  )

# -----------------------------------------------------------------------------
# Forest plot
# -----------------------------------------------------------------------------
plot_df <- bwqs_results
plot_df$Group <- factor(plot_df$Group, levels = rev(group_order))
plot_df$OR_CI_Label <- sprintf("%.2f (%.2f-%.2f)", plot_df$OR, plot_df$LCI, plot_df$UCI)
plot_df$Significance <- ifelse(
  plot_df$Model == "BWQS",
  ifelse(plot_df$LCI > 1 | plot_df$UCI < 1, "*", ""),
  ifelse(plot_df$Pvalue < 0.05, "*", "")
)

x_min <- min(c(0.4, min(plot_df$LCI) * 0.85), na.rm = TRUE)
x_max <- max(c(2.0, max(plot_df$UCI) * 1.2), na.rm = TRUE)

p <- ggplot(plot_df, aes(x = OR, y = Group)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 1.0) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI, color = Group),
                 height = 0.25, linewidth = 2.5) +
  geom_point(aes(fill = Group, color = Group),
             shape = 23, size = 6, stroke = 1.2) +
  geom_text(aes(x = UCI + (x_max - x_min) * 0.04, label = Significance),
            color = "red", size = 7, fontface = "bold",
            hjust = 0, vjust = 0.3, show.legend = FALSE) +
  geom_text(aes(x = x_max * 0.98, label = OR_CI_Label),
            hjust = 1, size = 4.2, color = "gray20") +
  scale_color_manual(values = user_colors, guide = "none") +
  scale_fill_manual(values = user_colors) +
  scale_x_continuous(
    name   = "Odds Ratio (95% CI)",
    limits = c(x_min, x_max * 1.18),
    expand = c(0, 0)
  ) +
  scale_y_discrete(name = "") +
  facet_grid(rows = vars(Model), scales = "free_y", space = "free_y",
             labeller = labeller(Model = c("BWQS" = "BWQS (Multi-pollutant)",
                                           "GLM"  = "GLM (Single-pollutant)"))) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "gray85", linewidth = 0.3),
    axis.text.y        = element_text(size = 12, face = "bold"),
    axis.text.x        = element_text(size = 11),
    axis.title.x       = element_text(size = 13, face = "bold"),
    strip.text         = element_text(size = 13, face = "bold"),
    strip.background   = element_rect(fill = "gray95", color = "gray70"),
    legend.position    = "bottom",
    legend.title       = element_text(size = 12, face = "bold"),
    legend.text        = element_text(size = 11),
    plot.title         = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle      = element_text(size = 11, hjust = 0.5, color = "gray40"),
    plot.caption       = element_text(size = 10, hjust = 1, color = "gray50"),
    plot.margin        = margin(t = 15, r = 50, b = 10, l = 10)
  ) +
  guides(fill = guide_legend(
    title = "Pollutant Group",
    ncol  = 4,
    override.aes = list(shape = 23, size = 5, stroke = 1)
  )) +
  labs(
    title    = "Association Between Pollutant Mixtures and Hypertension",
    subtitle = "Bayesian Weighted Quantile Sum (BWQS) & Single-Pollutant GLM",
    caption  = "* Significant: P < 0.05 (GLM) or 95% CI excludes 1 (BWQS)"
  )

ggsave("Forest_Plot_BWQS.pdf", plot = p,
       width = 16, height = max(6, 0.45 * nrow(plot_df) + 3),
       device = "pdf", dpi = 300)
ggsave("Forest_Plot_BWQS.png", plot = p,
       width = 16, height = max(6, 0.45 * nrow(plot_df) + 3),
       device = "png", dpi = 300)
