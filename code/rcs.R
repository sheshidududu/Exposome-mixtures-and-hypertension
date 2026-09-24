# =============================================================================
# Restricted cubic spline (RCS) dose-response analysis
# Single-pollutant logistic models with 3-knot RCS, per pollutant
# =============================================================================

rm(list = ls())

library(ggplot2)
library(rms)
library(openxlsx)

# -----------------------------------------------------------------------------
# Input (load before running)
# -----------------------------------------------------------------------------
# data_raw <- read.xlsx(...)

pollutants <- c("SFM", "4AN6C2B", "SFSX", "TMTH", "SFTHI",
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
                "Sr", "Mo", "Cd", "Te", "Tl", "Pb")

covariates_cat <- c("Gender", "Exercisefrequency", "Eatinghabits", "Smoke", "Drink")
covariates_cont <- c("Leftfastingbloodglucose2mmol", "Age2", "BMI")
all_covariates <- c(covariates_cat, covariates_cont)

outcome <- "Hypertension"

all_vars <- c(outcome, pollutants, all_covariates)
missing_vars <- all_vars[!all_vars %in% colnames(data_raw)]
if (length(missing_vars) > 0) {
  warning("Variables not in data: ", paste(missing_vars, collapse = ", "))
  pollutants <- pollutants[pollutants %in% colnames(data_raw)]
  cat("Pollutants available:", length(pollutants), "\n")
}

data <- na.omit(data_raw[, c(outcome, pollutants, all_covariates)])

data$Hypertension <- as.factor(data$Hypertension)
for (var in covariates_cat) {
  if (var %in% colnames(data)) data[[var]] <- as.factor(data[[var]])
}
for (var in covariates_cont) {
  if (var %in% colnames(data)) data[[var]] <- as.numeric(data[[var]])
}
for (var in pollutants) {
  if (var %in% colnames(data)) data[[var]] <- as.numeric(data[[var]])
}

results <- data.frame(
  Pollutant = character(),
  P_overall = numeric(),
  P_nonlinear = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(pollutants)) {
  pollutant <- pollutants[i]
  cat("Processing:", pollutant, "(", i, "/", length(pollutants), ")\n")

  core_name <- gsub("^ln_|_1$", "", pollutant)

  data_sub <- na.omit(data[, c(outcome, pollutant, all_covariates)])

  if (nrow(data_sub) < 50) {
    warning(paste("Sample size insufficient (", nrow(data_sub), "); skipped:", pollutant))
    results <- rbind(results, data.frame(
      Pollutant = pollutant, P_overall = NA, P_nonlinear = NA
    ))
    next
  }

  colnames(data_sub)[colnames(data_sub) == pollutant] <- "exposure"

  dd <- datadist(data_sub)
  options(datadist = "dd")

  formula_str <- paste("Hypertension ~ rcs(exposure, 3) +",
                       paste(all_covariates, collapse = " + "))
  model <- lrm(as.formula(formula_str), data = data_sub, x = TRUE, y = TRUE)

  if (any(is.na(coef(model)))) {
    warning(paste("Model failed (NA coefficients); skipped:", pollutant))
    results <- rbind(results, data.frame(
      Pollutant = pollutant, P_overall = NA, P_nonlinear = NA
    ))
    next
  }

  an <- anova(model)
  rownames_an <- rownames(an)
  idx_overall <- grep("exposure", rownames_an, ignore.case = TRUE)
  idx_overall <- idx_overall[!grepl("Nonlinear", rownames_an[idx_overall], ignore.case = TRUE)]
  idx_nonlinear <- grep("Nonlinear", rownames_an, ignore.case = TRUE)

  p_overall <- if (length(idx_overall) == 1) an[idx_overall, "P"] else NA
  p_nonlinear <- if (length(idx_nonlinear) == 1) an[idx_nonlinear, "P"] else NA

  OR_temp <- Predict(model, exposure, fun = exp, type = "predictions",
                     ref.zero = TRUE, np = 200)
  idx_ref <- which.min(abs(OR_temp$yhat - 1))
  ref_value <- OR_temp$exposure[idx_ref]
  dd$limits$exposure[2] <- ref_value
  options(datadist = "dd")
  model <- update(model)

  OR_final <- Predict(model, exposure, fun = exp, type = "predictions",
                      ref.zero = TRUE, np = 200)
  plot_df <- data.frame(
    exposure = OR_final$exposure,
    OR = OR_final$yhat,
    lower = OR_final$lower,
    upper = OR_final$upper
  )

  p <- ggplot(plot_df, aes(x = exposure, y = OR)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "grey40") +
    geom_line(color = "grey20", size = 1.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "orange", size = 2) +
    labs(x = "Concentration", y = "OR") +
    theme_bw(base_size = 16) +
    theme(
      aspect.ratio = 1,
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 20),
      axis.text.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
      axis.line = element_line(color = "black"),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = Inf, y = Inf, label = core_name,
             hjust = 1.1, vjust = 1.5, size = 15, fontface = "bold", color = "black")

  print(p)

  results <- rbind(results, data.frame(
    Pollutant = pollutant,
    P_overall = p_overall,
    P_nonlinear = p_nonlinear
  ))

  cat("  P-overall =", round(p_overall, 4),
      " | P-nonlinear =", round(p_nonlinear, 4), "\n")

  rm(data_sub, dd, model, OR_temp, OR_final, plot_df, p)
  gc()
}

results$P_overall <- round(results$P_overall, 3)
results$P_nonlinear <- round(results$P_nonlinear, 3)

cat("\nAll pollutants processed.\n")
