# =============================================================================
# Sensitivity analysis: log-transformed exposures, Winsorized at 1%/99%,
# single-pollutant logistic regression with and without eGFR adjustment
# (CKD-EPI 2021); Table S21 + forest plot
# =============================================================================

rm(list = ls())

library(readxl)
library(dplyr)
library(writexl)
library(ggplot2)
library(broom)
library(cowplot)

# -----------------------------------------------------------------------------
# Input (load before running)
# -----------------------------------------------------------------------------
# data <- read_excel(...)

data$Hypertension <- as.numeric(data$Hypertension)
data$Serumcreatininevalue2umolL <- as.numeric(data$Serumcreatininevalue2umolL)
data$Age2 <- as.numeric(data$Age2)
data$BMI <- as.numeric(data$BMI)
data$Leftfastingbloodglucose2mmol <- as.numeric(data$Leftfastingbloodglucose2mmol)

# -----------------------------------------------------------------------------
# eGFR (CKD-EPI 2021, creatinine in umol/L)
# -----------------------------------------------------------------------------
calculate_egfr <- function(scr, age, sex, units = "umol") {
  if (units == "umol") scr_mgdl <- scr / 88.4 else scr_mgdl <- scr
  if (sex == 1) { k <- 0.9; a <- 142 } else { k <- 0.7; a <- 144 }
  alpha1 <- -0.299; alpha2 <- -1.209
  scr_ratio <- scr_mgdl / k
  term1 <- ifelse(scr_ratio <= 1, scr_ratio^alpha1, scr_ratio^alpha2)
  a * term1 * (0.993^age)
}
data$eGFR <- mapply(calculate_egfr, data$Serumcreatininevalue2umolL,
                    data$Age2, data$Gender, units = "umol")
cat("eGFR Mean (SD):", round(mean(data$eGFR, na.rm = TRUE), 2),
    "(", round(sd(data$eGFR, na.rm = TRUE), 2), ")\n")

# -----------------------------------------------------------------------------
# Categorical covariates
# -----------------------------------------------------------------------------
categorical_vars <- c("Gender", "Exercisefrequency", "Eatinghabits", "Smoke", "Drink")
for (v in categorical_vars) {
  if (v %in% names(data)) data[[v]] <- as.factor(data[[v]])
}

# -----------------------------------------------------------------------------
# Pollutants (70, log-transformed)
# -----------------------------------------------------------------------------
excel_names <- c(
  # Sulfonamides (5)
  "SFM", "4AN6C2B", "SFSX", "TMTH", "SFTHI",
  # Fluoroquinolones (4)
  "NOR", "ENOR", "CIP", "OFLX",
  # Tetracyclines (4)
  "TET", "OTC", "AM", "DOXY",
  # Neonicotinoids (15)
  "IMI", "NIT", "ACE", "THM", "CLO", "THI", "IPP", "CYC", "DIN", "IMID",
  "DM-THM", "DN-IMI", "5-OH-IMI", "DM-ACE", "IMI-olefin",
  # Organophosphates (5)
  "DZN", "CPS-O", "3,5,6-TCP", "DDVP", "MALT",
  # Carbamates (1)
  "3OH-CBF",
  # Pyrethroids (5)
  "3-PBA", "4-F-3-PBA", "2-MB-3CA", "trans-DCCA", "cis-DCCA",
  # Herbicides (4)
  "ATZ", "ATZ-DE-DSI-2HD", "CBDZ", "DMM",
  # Phenylpyrazoles (1)
  "FP",
  # Bisphenols (5)
  "BPF", "BPS", "BPC", "BPE", "BPP",
  # OHPAHs (8)
  "2-OHN", "1-OHN", "3-OHF", "2-OHF", "2-OHP", "1-OHP", "4-OHP", "1-OHPYR",
  # Metals (13)
  "V", "Fe", "Co", "Cu", "Zn", "As", "Se", "Sr", "Mo", "Cd", "Te", "Tl", "Pb"
)
stopifnot(all(excel_names %in% names(data)))
cat("Pollutant count:", length(excel_names), "\n")

# -----------------------------------------------------------------------------
# Group labels and colors
# -----------------------------------------------------------------------------
group_labels <- c(
  rep("Sulfonamides", 5),
  rep("Fluoroquinolones", 4),
  rep("Tetracyclines", 4),
  rep("Neonicotinoids", 15),
  rep("Organophosphates", 5),
  "Carbamates",
  rep("Pyrethroids", 5),
  rep("Herbicides", 4),
  "Phenylpyrazoles",
  rep("Bisphenols", 5),
  rep("OHPAHs", 8),
  rep("Metals", 13)
)
stopifnot(length(excel_names) == length(group_labels))

group_colors <- c(
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
group_levels <- names(group_colors)

# -----------------------------------------------------------------------------
# Winsorize at 1%/99% (data already log-transformed)
# -----------------------------------------------------------------------------
winsorize_1_99 <- function(x) {
  x <- as.numeric(as.character(x))
  q_low <- quantile(x, 0.01, na.rm = TRUE)
  q_high <- quantile(x, 0.99, na.rm = TRUE)
  x[x < q_low] <- q_low
  x[x > q_high] <- q_high
  return(x)
}

for (p in excel_names) {
  data[[p]] <- winsorize_1_99(data[[p]])
}

# -----------------------------------------------------------------------------
# Complete cases
# -----------------------------------------------------------------------------
complete_vars <- c("Hypertension", categorical_vars, "Age2", "BMI",
                   "Leftfastingbloodglucose2mmol", "eGFR", excel_names)
df <- data[complete.cases(data[, complete_vars]), ]
cat("Sample:", nrow(df), "  Case:", sum(df$Hypertension == 1),
    "  Control:", sum(df$Hypertension == 0), "\n")

# -----------------------------------------------------------------------------
# Two models: without eGFR vs with eGFR
# -----------------------------------------------------------------------------
cov_no <- c("Gender", "Age2", "BMI", "Exercisefrequency", "Eatinghabits",
            "Smoke", "Drink", "Leftfastingbloodglucose2mmol")
cov_with <- c(cov_no, "eGFR")

results_list <- list()
n_fail <- 0

for (i in seq_along(excel_names)) {
  p <- excel_names[i]
  grp <- group_labels[i]

  f1 <- paste0("Hypertension ~ `", p, "` + ", paste(cov_no, collapse = "+"))
  f2 <- paste0("Hypertension ~ `", p, "` + ", paste(cov_with, collapse = "+"))

  mod1 <- try(glm(as.formula(f1), data = df, family = binomial), silent = TRUE)
  mod2 <- try(glm(as.formula(f2), data = df, family = binomial), silent = TRUE)

  if (inherits(mod1, "try-error") || inherits(mod2, "try-error")) {
    n_fail <- n_fail + 1
    cat(sprintf("  x %-22s | model fitting failed\n", p))
    next
  }

  term_p <- paste0("`", p, "`")
  r1 <- tidy(mod1, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == p | term == term_p)
  r2 <- tidy(mod2, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == p | term == term_p)

  if (nrow(r1) == 0 || nrow(r2) == 0) {
    n_fail <- n_fail + 1
    cat(sprintf("  x %-22s | tidy failed\n", p))
    next
  }

  or1 <- r1$estimate; lo1 <- r1$conf.low; hi1 <- r1$conf.high; pv1 <- r1$p.value
  or2 <- r2$estimate; lo2 <- r2$conf.low; hi2 <- r2$conf.high; pv2 <- r2$p.value
  delta <- (or2 - or1) / or1 * 100

  results_list[[length(results_list) + 1]] <- data.frame(
    Pollutant = p, Group = grp,
    OR_no = or1, CI_low_no = lo1, CI_up_no = hi1, P_no = pv1,
    OR_with = or2, CI_low_with = lo2, CI_up_with = hi2, P_with = pv2,
    OR_change_pct = delta,
    stringsAsFactors = FALSE
  )
  cat(sprintf("  ok %-22s | OR_no: %.4f | OR_with: %.4f | delta: %+.3f%%\n",
              p, or1, or2, delta))
}
cat(sprintf("Success: %d / %d   Failed: %d\n",
            length(results_list), length(excel_names), n_fail))
stopifnot(length(results_list) > 0)

res_df <- bind_rows(results_list)
res_df$Group <- factor(res_df$Group, levels = group_levels)

# -----------------------------------------------------------------------------
# Table S21
# -----------------------------------------------------------------------------
table_s21 <- res_df %>%
  mutate(
    `OR (without eGFR)` = sprintf("%.3f (%.3f, %.3f)", OR_no, CI_low_no, CI_up_no),
    `P (without eGFR)`  = ifelse(P_no < 0.001, "<0.001", sprintf("%.3f", P_no)),
    `OR (with eGFR)`    = sprintf("%.3f (%.3f, %.3f)", OR_with, CI_low_with, CI_up_with),
    `P (with eGFR)`     = ifelse(P_with < 0.001, "<0.001", sprintf("%.3f", P_with)),
    `OR change (%)`     = sprintf("%.2f", OR_change_pct)
  ) %>%
  select(Pollutant, Group,
         `OR (without eGFR)`, `P (without eGFR)`,
         `OR (with eGFR)`, `P (with eGFR)`,
         `OR change (%)`)

# -----------------------------------------------------------------------------
# Forest plot (with-eGFR model)
# -----------------------------------------------------------------------------
plot_data <- res_df %>%
  rename(OR = OR_with, Lower = CI_low_with, Upper = CI_up_with, P_val = P_with) %>%
  mutate(
    FDR        = p.adjust(P_val, method = "BH"),
    Group_Hex  = group_colors[as.character(Group)],
    Point_Fill = ifelse(P_val < 0.05, Group_Hex, "white"),
    Star_Label = ifelse(FDR < 0.05, "*", NA)
  )

create_main_plot <- function(data_all) {
  data_all$Label <- factor(data_all$Pollutant, levels = rev(excel_names))
  data_all$Facet <- "Sensitivity: Winsorized + eGFR Adjusted"

  min_val <- min(data_all$Lower, na.rm = TRUE)
  max_val <- max(data_all$Upper, na.rm = TRUE)
  x_min <- min_val * 0.9; x_max <- max_val * 1.15
  rect_w <- (x_max - x_min) * 0.04

  ggplot(data_all, aes(y = Label)) +
    geom_rect(aes(xmin = x_min - rect_w, xmax = x_min,
                  ymin = as.numeric(Label) - 0.5, ymax = as.numeric(Label) + 0.5,
                  fill = Group), color = "black", linewidth = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_errorbar(aes(xmin = Lower, xmax = Upper, color = Group),
                  linewidth = 1.5, width = 0.45) +
    geom_point(aes(x = OR, color = Group, fill = I(Point_Fill)),
               shape = 23, size = 4.5, stroke = 0.7) +
    geom_text(aes(x = Upper, label = Star_Label),
              hjust = -0.5, vjust = 0.7, size = 5.5, fontface = "bold", na.rm = TRUE) +
    facet_grid(. ~ Facet) +
    scale_x_continuous(limits = c(x_min - rect_w, x_max), trans = "log10",
                       expand = c(0, 0), breaks = c(1)) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    labs(x = "OR (95% CI)", y = NULL) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.y  = element_text(size = 9, color = "black", face = "bold"),
      axis.text.x  = element_text(size = 11),
      axis.title.x = element_text(size = 12, face = "bold"),
      panel.grid   = element_blank(),
      strip.background = element_rect(fill = "#2A6A72", color = "black"),
      strip.text   = element_text(color = "white", face = "bold", size = 12),
      plot.margin  = margin(10, 20, 5, 10)
    )
}

get_legend <- function(data_all) {
  p1 <- ggplot(data_all, aes(x = OR, y = Pollutant, fill = Group)) +
    geom_point(shape = 22, size = 5, color = "transparent") +
    scale_fill_manual(values = group_colors) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 8.5, face = "bold")) +
    guides(fill = guide_legend(nrow = 3))
  leg1 <- cowplot::get_legend(p1)

  sig_df <- data.frame(
    Label = factor(c("Significant (P<0.05)", "Non-significant")),
    x = 1:2, y = c(1, 1)
  )
  p2 <- ggplot(sig_df, aes(x, y, fill = Label)) +
    geom_point(shape = 23, size = 4, stroke = 1, color = "gray40") +
    scale_fill_manual(values = c("Significant (P<0.05)" = "gray40",
                                 "Non-significant" = "white")) +
    theme(legend.position = "bottom", legend.key = element_blank(),
          legend.text = element_text(size = 9, face = "bold"))
  leg2 <- cowplot::get_legend(p2)

  leg3 <- ggdraw() + draw_label("* FDR < 0.05", fontface = "bold", size = 9, x = 0.5, y = 0.6)
  plot_grid(leg1, leg2, leg3, ncol = 1, rel_heights = c(1.2, 0.4, 0.2))
}

plot_data_sorted <- plot_data[match(excel_names, plot_data$Pollutant), ]
p_main <- create_main_plot(plot_data_sorted)
p_leg <- get_legend(plot_data_sorted)
p_final <- plot_grid(p_main, p_leg, ncol = 1, rel_heights = c(1, 0.18))

ggsave("Fig_Sensitivity_Forest.pdf", p_final, width = 10, height = 16)
ggsave("Fig_Sensitivity_Forest.png", p_final, width = 10, height = 16, dpi = 300)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n--- OR change (%) ---\n")
cat(sprintf("  Mean +/- SD: %.3f +/- %.3f\n",
            mean(res_df$OR_change_pct), sd(res_df$OR_change_pct)))
cat(sprintf("  Range: [%.3f, %.3f]\n",
            min(res_df$OR_change_pct), max(res_df$OR_change_pct)))
cat(sprintf("  |delta|>1%%: %d / %d\n",
            sum(abs(res_df$OR_change_pct) > 1), nrow(res_df)))

cat("\n--- FDR<0.05 (with eGFR) ---\n")
sig <- plot_data %>% filter(FDR < 0.05) %>% arrange(FDR)
if (nrow(sig) > 0) {
  print(as.data.frame(sig %>% select(Pollutant, Group, OR, Lower, Upper, P_val, FDR)))
}
