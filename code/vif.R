# =============================================================================
# Variance inflation factor (VIF) collinearity diagnostics
# Full-model correlation-matrix method: pre-merge vs post-merge
# Merged pairs: trans-DCCA + cis-DCCA -> DCCA; TET + DOXY -> TETDOXY
# =============================================================================

rm(list = ls())

library(openxlsx)
library(dplyr)

# -----------------------------------------------------------------------------
# Input (load before running)
# -----------------------------------------------------------------------------
# data <- read.xlsx(...)

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
# Merge highly collinear variables
# -----------------------------------------------------------------------------
high_corr_pairs <- list(
  DCCA    = c("trans-DCCA", "cis-DCCA"),
  TETDOXY = c("TET", "DOXY")
)

for (pair_name in names(high_corr_pairs)) {
  vars <- high_corr_pairs[[pair_name]]
  if (all(vars %in% names(data))) {
    data[[pair_name]] <- rowMeans(data[, vars], na.rm = TRUE)
    cat(sprintf("Merged variable: %s = mean(%s)\n", pair_name, paste(vars, collapse = " + ")))
  }
}

# -----------------------------------------------------------------------------
# Pre-merge / post-merge full-model predictor lists
# -----------------------------------------------------------------------------
all_pollutants_pre <- unique(unlist(pollutant_groups))
all_pollutants_pre <- intersect(all_pollutants_pre, names(data))

all_pollutants_post <- setdiff(all_pollutants_pre,
                               unlist(high_corr_pairs, use.names = FALSE))
all_pollutants_post <- c(all_pollutants_post, names(high_corr_pairs))
all_pollutants_post <- intersect(all_pollutants_post, names(data))

cat(sprintf("\nPre-merge full model: %d pollutants + %d covariates\n",
            length(all_pollutants_pre), length(covariates)))
cat(sprintf("Post-merge full model: %d pollutants + %d covariates\n",
            length(all_pollutants_post), length(covariates)))

# -----------------------------------------------------------------------------
# Full-model VIF (correlation matrix method): VIF_j = diag(R^-1)_j
# -----------------------------------------------------------------------------
compute_full_vif <- function(pollutant_list, df, label) {
  cat(sprintf("\n--- %s ---\n", label))

  num_covariates <- c("Age2", "BMI", "Leftfastingbloodglucose2mmol")

  vars_all <- c(pollutant_list, covariates)
  vars_exist <- intersect(vars_all, names(df))
  df_sub <- na.omit(df[, vars_exist])

  safe <- function(x) sapply(x, function(v) {
    if (grepl("^[a-zA-Z][a-zA-Z0-9._]*$", v)) v else paste0("`", v, "`")
  }, USE.NAMES = FALSE)
  form <- as.formula(paste("~", paste(safe(vars_exist), collapse = " + ")))

  X <- model.matrix(form, data = df_sub)
  X <- X[, -1, drop = FALSE]

  cat(sprintf("  Complete cases: %d, design matrix columns: %d\n", nrow(X), ncol(X)))

  X_scaled <- scale(X, center = TRUE, scale = TRUE)
  R <- cor(X_scaled, use = "complete.obs")

  eig <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(eig)
  cat(sprintf("  Minimum eigenvalue of correlation matrix: %.4e\n", min_eig))

  is_singular <- min_eig < 1e-10

  if (!is_singular) {
    R_inv <- solve(R)
    vif_all <- diag(R_inv)
  } else {
    # Truncated-eigenvalue pseudo-inverse
    eig_full <- eigen(R, symmetric = TRUE)
    pos <- which(eig_full$values > 1e-10)
    cat(sprintf("  Matrix nearly singular, effective rank %d/%d, pseudo-inverse used\n",
                length(pos), ncol(R)))
    R_inv_approx <- eig_full$vectors[, pos] %*%
      diag(1 / eig_full$values[pos], nrow = length(pos)) %*%
      t(eig_full$vectors[, pos])
    vif_all <- diag(R_inv_approx)
  }
  names(vif_all) <- colnames(R)

  result <- data.frame(Variable = pollutant_list, VIF = NA_real_,
                       stringsAsFactors = FALSE)
  for (i in seq_len(nrow(result))) {
    v <- result$Variable[i]
    if (v %in% names(vif_all)) {
      result$VIF[i] <- round(vif_all[v], 4)
    }
  }
  return(result)
}

# -----------------------------------------------------------------------------
# Compute pre- and post-merge VIF
# -----------------------------------------------------------------------------
vif_pre <- compute_full_vif(all_pollutants_pre, data, "Pre-Merge")
vif_post <- compute_full_vif(all_pollutants_post, data, "Post-Merge")

# -----------------------------------------------------------------------------
# Combine into one table
# -----------------------------------------------------------------------------
pollutant_to_group <- data.frame(Variable = character(), Group = character(),
                                 stringsAsFactors = FALSE)
for (grp in names(pollutant_groups)) {
  for (vv in pollutant_groups[[grp]]) {
    pollutant_to_group <- rbind(pollutant_to_group,
                                data.frame(Variable = vv, Group = grp,
                                           stringsAsFactors = FALSE))
  }
}
for (pair_name in names(high_corr_pairs)) {
  grp_found <- NA_character_
  for (grp in names(pollutant_groups)) {
    if (any(high_corr_pairs[[pair_name]] %in% pollutant_groups[[grp]])) {
      grp_found <- grp
      break
    }
  }
  pollutant_to_group <- rbind(pollutant_to_group,
                              data.frame(Variable = pair_name, Group = grp_found,
                                         stringsAsFactors = FALSE))
}

colnames(vif_pre) <- c("Variable", "VIF_PreMerge")
colnames(vif_post) <- c("Variable", "VIF_PostMerge")

combined <- merge(vif_pre, vif_post, by = "Variable", all = TRUE)
combined <- merge(combined, pollutant_to_group, by = "Variable", all.x = TRUE)
combined$Group <- factor(combined$Group, levels = group_order)

combined$Merge_Status <- sapply(combined$Variable, function(v) {
  if (v %in% names(high_corr_pairs)) {
    paste0("<- merged: ", paste(high_corr_pairs[[v]], collapse = " + "))
  } else if (v %in% unlist(high_corr_pairs, use.names = FALSE)) {
    "replaced by merged variable"
  } else {
    "-"
  }
})

combined$VIF_Delta <- round(combined$VIF_PostMerge - combined$VIF_PreMerge, 4)

judge <- function(v) {
  ifelse(is.na(v), "-",
         ifelse(v < 5, "Good (<5)",
                ifelse(v < 10, "Moderate (5-10)",
                       "Severe (>10)")))
}
combined$PreMerge_Judge <- judge(combined$VIF_PreMerge)
combined$PostMerge_Judge <- judge(combined$VIF_PostMerge)

combined <- combined %>% arrange(Group, Variable)

colnames(combined) <- c(
  "Variable", "VIF (Pre-Merge)", "VIF (Post-Merge)",
  "Chemical Group", "Merge Status", "VIF Delta (Post-Pre)",
  "Pre-Merge Judgment", "Post-Merge Judgment"
)

# -----------------------------------------------------------------------------
# Key results
# -----------------------------------------------------------------------------
cat(sprintf("\nPre-merge (%d pollutants):\n", nrow(vif_pre)))
cat(sprintf("  VIF > 10 (severe): %d\n", sum(vif_pre$VIF > 10, na.rm = TRUE)))
cat(sprintf("  VIF 5-10 (moderate): %d\n", sum(vif_pre$VIF >= 5 & vif_pre$VIF <= 10, na.rm = TRUE)))
cat(sprintf("  VIF < 5 (good): %d\n", sum(vif_pre$VIF < 5, na.rm = TRUE)))
cat(sprintf("  VIF = NA: %d\n", sum(is.na(vif_pre$VIF))))

cat("\nPost-merge (%d pollutants):\n")
cat(sprintf("  VIF > 10 (severe): %d\n", sum(vif_post$VIF > 10, na.rm = TRUE)))
cat(sprintf("  VIF 5-10 (moderate): %d\n", sum(vif_post$VIF >= 5 & vif_post$VIF <= 10, na.rm = TRUE)))
cat(sprintf("  VIF < 5 (good): %d\n", sum(vif_post$VIF < 5, na.rm = TRUE)))
cat(sprintf("  VIF = NA: %d\n", sum(is.na(vif_post$VIF))))

cat("\n--- VIF change for merged variables ---\n")
for (pair_name in names(high_corr_pairs)) {
  old_vars <- high_corr_pairs[[pair_name]]
  cat(sprintf("\n> %s = mean(%s):\n", pair_name, paste(old_vars, collapse = ", ")))
  for (ov in old_vars) {
    ovif <- vif_pre$VIF[vif_pre$Variable == ov]
    cat(sprintf("    %-20s pre-merge VIF = %s\n", ov,
                if (length(ovif) > 0 && !is.na(ovif)) sprintf("%.4f", ovif) else "-"))
  }
  nvif <- vif_post$VIF[vif_post$Variable == pair_name]
  cat(sprintf("    -> %-18s post-merge VIF = %s\n", pair_name,
              if (length(nvif) > 0 && !is.na(nvif)) sprintf("%.4f", nvif) else "-"))
}

cat("\n--- Variables with pre-merge VIF >= 5 ---\n")
high_vif <- combined[combined$`VIF (Pre-Merge)` >= 5 & !is.na(combined$`VIF (Pre-Merge)`), ]
if (nrow(high_vif) > 0) {
  for (i in seq_len(nrow(high_vif))) {
    cat(sprintf("  %-20s  VIF=%.2f  [%s]\n",
                high_vif$Variable[i], high_vif$`VIF (Pre-Merge)`[i],
                high_vif$`Chemical Group`[i]))
  }
} else {
  cat("  (none)\n")
}
