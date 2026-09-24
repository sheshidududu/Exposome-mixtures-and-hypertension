# =============================================================================
# Quartile logistic regression (unadjusted / adjusted) with FDR correction,
# plus continuous-variable logistic regression with FDR correction
# Table output only, no figures
# =============================================================================

rm(list = ls())

library(readxl)
library(dplyr)
library(writexl)

# -----------------------------------------------------------------------------
# Input (load before running)
# -----------------------------------------------------------------------------
# data <- read_excel(...)

data$Hypertension <- as.numeric(data$Hypertension)
cat("Outcome distribution:\n")
print(table(data$Hypertension, useNA = "ifany"))

# -----------------------------------------------------------------------------
# Pollutants (70)
# -----------------------------------------------------------------------------
pollutants <- c(
  "SFM", "4AN6C2B", "SFSX", "TMTH", "SFTHI", "NOR", "ENOR", "CIP", "OFLX", "TET",
  "OTC", "AM", "DOXY", "IMI", "NIT", "ACE", "THM", "CLO", "THI", "IPP", "CYC",
  "DIN", "IMID", "DM-THM", "DN-IMI", "5-OH-IMI", "DM-ACE", "FP", "DMM", "CBDZ",
  "CPS-O", "DZN", "DDVP", "MALT", "3OH-CBF", "ATZ", "ATZ-DE-DSI-2HD",
  "BPF", "BPS", "BPC", "BPE", "BPP",
  "2-OHN", "1-OHN", "3-OHF", "2-OHF", "2-OHP", "1-OHP", "4-OHP", "1-OHPYR",
  "3-PBA", "4-F-3-PBA", "2-MB-3CA", "trans-DCCA", "cis-DCCA", "IMI-olefin", "3,5,6-TCP",
  "V", "Fe", "Co", "Cu", "Zn", "As", "Se", "Sr", "Mo", "Cd", "Te", "Tl", "Pb"
)

pollutants <- intersect(pollutants, names(data))
cat("Pollutants available:", length(pollutants), "\n")
print(pollutants)

# -----------------------------------------------------------------------------
# Covariates
# -----------------------------------------------------------------------------
categorical_vars <- c("Gender", "Exercisefrequency", "Eatinghabits", "Smoke", "Drink")
for (v in categorical_vars) {
  if (v %in% names(data)) data[[v]] <- as.factor(data[[v]])
}
data$Age2 <- as.numeric(data$Age2)
data$BMI <- as.numeric(data$BMI)
data$Leftfastingbloodglucose2mmol <- as.numeric(data$Leftfastingbloodglucose2mmol)

for (p in pollutants) {
  data[[p]] <- as.numeric(as.character(data[[p]]))
}

# -----------------------------------------------------------------------------
# Complete cases
# -----------------------------------------------------------------------------
complete_vars <- c("Hypertension", categorical_vars, "Age2", "BMI",
                   "Leftfastingbloodglucose2mmol", pollutants)
data_comp <- data[complete.cases(data[, complete_vars]), ]
cat("Complete sample:", nrow(data_comp), "\n")
if (nrow(data_comp) < 10) stop("Sample size too small")

# -----------------------------------------------------------------------------
# Quartile models
# -----------------------------------------------------------------------------
unadj_results <- list()
adj_results <- list()

for (p in pollutants) {
  cat("\nProcessing:", p, "\n")

  x <- data_comp[[p]]
  if (length(unique(x)) < 4) {
    cat("  -> fewer than 4 unique values; skipped\n")
    next
  }

  quants <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  if (length(unique(quants)) < 4) {
    cat("  -> fewer than 4 distinct quartile breaks; skipped\n")
    next
  }

  grp <- cut(x, breaks = quants, include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))
  if (any(table(grp) == 0)) {
    cat("  -> empty group; skipped\n")
    next
  }
  grp <- relevel(grp, ref = "Q1")

  temp <- data.frame(
    y = data_comp$Hypertension,
    grp = grp,
    Gender = data_comp$Gender,
    Age2 = data_comp$Age2,
    BMI = data_comp$BMI,
    Exercisefrequency = data_comp$Exercisefrequency,
    Eatinghabits = data_comp$Eatinghabits,
    Smoke = data_comp$Smoke,
    Drink = data_comp$Drink,
    Leftfastingbloodglucose2mmol = data_comp$Leftfastingbloodglucose2mmol
  )

  # Unadjusted quartile ORs
  mod1 <- try(glm(y ~ grp, data = temp, family = binomial), silent = TRUE)
  if (!inherits(mod1, "try-error")) {
    coef1 <- summary(mod1)$coefficients
    for (lev in c("Q2", "Q3", "Q4")) {
      rn <- paste0("grp", lev)
      if (rn %in% rownames(coef1)) {
        beta <- coef1[rn, "Estimate"]
        pval <- coef1[rn, "Pr(>|z|)"]
        or <- exp(beta)
        ci <- try(exp(confint.default(mod1))[rn, ], silent = TRUE)
        if (!inherits(ci, "try-error") && !any(is.na(ci))) {
          unadj_results[[length(unadj_results) + 1]] <- data.frame(
            Pollutant = p, Group = lev,
            OR = or, CI_low = ci[1], CI_up = ci[2],
            P = pval, stringsAsFactors = FALSE
          )
        }
      }
    }

    # Unadjusted trend test
    temp$grp_num <- as.numeric(temp$grp)
    mod_trend <- try(glm(y ~ grp_num, data = temp, family = binomial), silent = TRUE)
    if (!inherits(mod_trend, "try-error")) {
      coef_t <- summary(mod_trend)$coefficients
      if (nrow(coef_t) >= 2) {
        pval_t <- coef_t[2, "Pr(>|z|)"]
        unadj_results[[length(unadj_results) + 1]] <- data.frame(
          Pollutant = p, Group = "P for trend",
          OR = NA, CI_low = NA, CI_up = NA,
          P = pval_t, stringsAsFactors = FALSE
        )
      }
    }
  } else {
    cat("  -> unadjusted model failed\n")
  }

  # Adjusted quartile ORs
  mod2 <- try(glm(y ~ grp + Gender + Age2 + BMI + Exercisefrequency + Eatinghabits +
                    Smoke + Drink + Leftfastingbloodglucose2mmol,
                  data = temp, family = binomial), silent = TRUE)
  if (!inherits(mod2, "try-error")) {
    coef2 <- summary(mod2)$coefficients
    for (lev in c("Q2", "Q3", "Q4")) {
      rn <- paste0("grp", lev)
      if (rn %in% rownames(coef2)) {
        beta <- coef2[rn, "Estimate"]
        pval <- coef2[rn, "Pr(>|z|)"]
        or <- exp(beta)
        ci <- try(exp(confint.default(mod2))[rn, ], silent = TRUE)
        if (!inherits(ci, "try-error") && !any(is.na(ci))) {
          adj_results[[length(adj_results) + 1]] <- data.frame(
            Pollutant = p, Group = lev,
            OR = or, CI_low = ci[1], CI_up = ci[2],
            P = pval, stringsAsFactors = FALSE
          )
        }
      }
    }

    # Adjusted trend test
    temp$grp_num <- as.numeric(temp$grp)
    mod_trend2 <- try(glm(y ~ grp_num + Gender + Age2 + BMI + Exercisefrequency +
                            Eatinghabits + Smoke + Drink + Leftfastingbloodglucose2mmol,
                          data = temp, family = binomial), silent = TRUE)
    if (!inherits(mod_trend2, "try-error")) {
      coef_t2 <- summary(mod_trend2)$coefficients
      if (nrow(coef_t2) >= 2) {
        pval_t2 <- coef_t2[2, "Pr(>|z|)"]
        adj_results[[length(adj_results) + 1]] <- data.frame(
          Pollutant = p, Group = "P for trend",
          OR = NA, CI_low = NA, CI_up = NA,
          P = pval_t2, stringsAsFactors = FALSE
        )
      }
    }
  } else {
    cat("  -> adjusted model failed\n")
  }
}

if (length(unadj_results) == 0) stop("No unadjusted results")
if (length(adj_results) == 0) stop("No adjusted results")

df_un <- do.call(rbind, unadj_results)
df_ad <- do.call(rbind, adj_results)

# -----------------------------------------------------------------------------
# FDR correction and formatting
# -----------------------------------------------------------------------------
df_un$FDR <- p.adjust(df_un$P, method = "fdr")
df_ad$FDR <- p.adjust(df_ad$P, method = "fdr")

df_un$OR_CI <- ifelse(is.na(df_un$OR), "", sprintf("%.3f (%.3f-%.3f)", df_un$OR, df_un$CI_low, df_un$CI_up))
df_ad$OR_CI <- ifelse(is.na(df_ad$OR), "", sprintf("%.3f (%.3f-%.3f)", df_ad$OR, df_ad$CI_low, df_ad$CI_up))

final_un <- df_un[, c("Pollutant", "Group", "OR_CI", "P", "FDR")]
final_ad <- df_ad[, c("Pollutant", "Group", "OR_CI", "P", "FDR")]

final_un$P <- round(final_un$P, 4)
final_un$FDR <- round(final_un$FDR, 4)
final_ad$P <- round(final_ad$P, 4)
final_ad$FDR <- round(final_ad$FDR, 4)

cat("\n\n================== Quartile unadjusted results ==================\n")
print(as.data.frame(final_un))

cat("\n\n================== Quartile adjusted results ==================\n")
print(as.data.frame(final_ad))

# -----------------------------------------------------------------------------
# Continuous-variable logistic regression (raw pollutant concentration)
# -----------------------------------------------------------------------------
cont_unadj_results <- list()
cont_adj_results <- list()

for (p in pollutants) {
  cat("\nContinuous variable:", p, "\n")

  temp <- data.frame(
    y = data_comp$Hypertension,
    x = data_comp[[p]],
    Gender = data_comp$Gender,
    Age2 = data_comp$Age2,
    BMI = data_comp$BMI,
    Exercisefrequency = data_comp$Exercisefrequency,
    Eatinghabits = data_comp$Eatinghabits,
    Smoke = data_comp$Smoke,
    Drink = data_comp$Drink,
    Leftfastingbloodglucose2mmol = data_comp$Leftfastingbloodglucose2mmol
  )

  # Unadjusted continuous model
  mod_cont1 <- try(glm(y ~ x, data = temp, family = binomial), silent = TRUE)
  if (!inherits(mod_cont1, "try-error")) {
    coef_c1 <- summary(mod_cont1)$coefficients
    if ("x" %in% rownames(coef_c1)) {
      beta <- coef_c1["x", "Estimate"]
      pval <- coef_c1["x", "Pr(>|z|)"]
      or <- exp(beta)
      ci <- try(exp(confint.default(mod_cont1))["x", ], silent = TRUE)
      if (!inherits(ci, "try-error") && !any(is.na(ci))) {
        cont_unadj_results[[length(cont_unadj_results) + 1]] <- data.frame(
          Pollutant = p,
          OR = or, CI_low = ci[1], CI_up = ci[2],
          P = pval, stringsAsFactors = FALSE
        )
      }
    }
  }

  # Adjusted continuous model
  mod_cont2 <- try(glm(y ~ x + Gender + Age2 + BMI + Exercisefrequency +
                         Eatinghabits + Smoke + Drink + Leftfastingbloodglucose2mmol,
                       data = temp, family = binomial), silent = TRUE)
  if (!inherits(mod_cont2, "try-error")) {
    coef_c2 <- summary(mod_cont2)$coefficients
    if ("x" %in% rownames(coef_c2)) {
      beta <- coef_c2["x", "Estimate"]
      pval <- coef_c2["x", "Pr(>|z|)"]
      or <- exp(beta)
      ci <- try(exp(confint.default(mod_cont2))["x", ], silent = TRUE)
      if (!inherits(ci, "try-error") && !any(is.na(ci))) {
        cont_adj_results[[length(cont_adj_results) + 1]] <- data.frame(
          Pollutant = p,
          OR = or, CI_low = ci[1], CI_up = ci[2],
          P = pval, stringsAsFactors = FALSE
        )
      }
    }
  }
}

df_cont_un <- do.call(rbind, cont_unadj_results)
df_cont_ad <- do.call(rbind, cont_adj_results)

df_cont_un$FDR <- p.adjust(df_cont_un$P, method = "fdr")
df_cont_ad$FDR <- p.adjust(df_cont_ad$P, method = "fdr")

df_cont_un$OR_CI <- sprintf("%.3f (%.3f-%.3f)", df_cont_un$OR, df_cont_un$CI_low, df_cont_un$CI_up)
df_cont_ad$OR_CI <- sprintf("%.3f (%.3f-%.3f)", df_cont_ad$OR, df_cont_ad$CI_low, df_cont_ad$CI_up)

final_cont_un <- df_cont_un[, c("Pollutant", "OR_CI", "P", "FDR")]
final_cont_ad <- df_cont_ad[, c("Pollutant", "OR_CI", "P", "FDR")]

final_cont_un$P <- round(final_cont_un$P, 4)
final_cont_un$FDR <- round(final_cont_un$FDR, 4)
final_cont_ad$P <- round(final_cont_ad$P, 4)
final_cont_ad$FDR <- round(final_cont_ad$FDR, 4)

cat("\n\n================== Continuous unadjusted results ==================\n")
print(as.data.frame(final_cont_un))

cat("\n\n================== Continuous adjusted results ==================\n")
print(as.data.frame(final_cont_ad))
