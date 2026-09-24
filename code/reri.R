# =============================================================================
# Pairwise interaction analysis: RERI, AP, S (delta method + FDR correction)
# =============================================================================

if (!require("openxlsx")) install.packages("openxlsx")
if (!require("dplyr")) install.packages("dplyr")
if (!require("pbapply")) install.packages("pbapply")

library(openxlsx)
library(dplyr)
library(pbapply)

# -----------------------------------------------------------------------------
# Input (load before running)
# -----------------------------------------------------------------------------
# data_raw <- read.xlsx(...)
# orig_colnames <- colnames(data_raw)
# valid_colnames <- make.names(orig_colnames, unique = TRUE)
# colnames(data_raw) <- valid_colnames
# name_map <- setNames(orig_colnames, valid_colnames)

# -----------------------------------------------------------------------------
# Pollutants (70) and covariates
# -----------------------------------------------------------------------------
pollutants_orig <- c("SFM", "4AN6C2B", "SFSX", "TMTH", "SFTHI", "NOR", "ENOR", "CIP", "OFLX",
                     "TET", "OTC", "AM", "DOXY", "IMI", "NIT", "ACE", "THM", "CLO", "THI",
                     "IPP", "CYC", "DIN", "IMID", "DM-THM", "DN-IMI", "5-OH-IMI", "DM-ACE",
                     "FP", "DMM", "CBDZ", "CPS-O", "DZN", "DDVP", "MALT", "3OH-CBF", "ATZ",
                     "ATZ-DE-DSI-2HD", "BPF",
                     "BPS", "BPC", "BPE", "BPP", "2-OHN", "1-OHN", "3-OHF", "2-OHF", "2-OHP",
                     "1-OHP", "4-OHP", "1-OHPYR", "3-PBA", "4-F-3-PBA", "2-MB-3CA",
                     "trans-DCCA", "cis-DCCA", "IMI-olefin", "3,5,6-TCP", "V", "Fe", "Co",
                     "Cu", "Zn", "As", "Se", "Sr", "Mo", "Cd", "Te", "Tl", "Pb")

pollutants_valid <- make.names(pollutants_orig)
pollutants_valid <- intersect(pollutants_valid, colnames(data_raw))
pollutants_orig <- name_map[pollutants_valid]
names(pollutants_valid) <- pollutants_orig

cat("Pollutants used:", length(pollutants_valid), "\n")

covariates_orig <- c("Gender", "Age2", "BMI", "Exercisefrequency",
                     "Eatinghabits", "Smoke", "Drink", "Leftfastingbloodglucose2mmol")
covariates_valid <- make.names(covariates_orig)
covariates_valid <- intersect(covariates_valid, colnames(data_raw))
covariates_orig <- name_map[covariates_valid]

for (cov in covariates_valid) {
  if (is.character(data_raw[[cov]]) || is.factor(data_raw[[cov]])) {
    data_raw[[cov]] <- as.factor(data_raw[[cov]])
  }
}

# -----------------------------------------------------------------------------
# Dichotomise pollutants at the median
# -----------------------------------------------------------------------------
data_bin <- data_raw
for (p in pollutants_valid) {
  med <- median(data_raw[[p]], na.rm = TRUE)
  data_bin[[p]] <- ifelse(data_raw[[p]] > med, 1, 0)
}

if (!is.numeric(data_bin$Hypertension)) {
  data_bin$Hypertension <- as.numeric(as.factor(data_bin$Hypertension)) - 1
}

# -----------------------------------------------------------------------------
# Pairwise interaction (RERI / AP / S)
# -----------------------------------------------------------------------------
calc_interaction <- function(p1_valid, p2_valid, data, outcome, covs_valid) {
  form_str <- paste(outcome, "~", p1_valid, "*", p2_valid, "+",
                    paste(covs_valid, collapse = "+"))
  form <- as.formula(form_str)

  fit <- tryCatch(
    glm(form, data = data, family = binomial),
    error = function(e) return(NULL)
  )
  if (is.null(fit)) return(NULL)

  coefs <- coef(fit)
  vcov_mat <- vcov(fit)

  idx1 <- match(p1_valid, names(coefs))
  idx2 <- match(p2_valid, names(coefs))
  idx_int <- match(paste0(p1_valid, ":", p2_valid), names(coefs))
  if (is.na(idx_int)) idx_int <- match(paste0(p2_valid, ":", p1_valid), names(coefs))

  if (any(is.na(c(idx1, idx2, idx_int)))) return(NULL)

  b1 <- coefs[idx1]
  b2 <- coefs[idx2]
  b3 <- coefs[idx_int]

  OR10 <- exp(b1)
  OR01 <- exp(b2)
  OR11 <- exp(b1 + b2 + b3)

  RERI <- OR11 - OR10 - OR01 + 1
  AP <- RERI / OR11
  denom <- (OR10 - 1) + (OR01 - 1)
  S <- ifelse(denom == 0, NA, (OR11 - 1) / denom)

  # Delta method standard error
  G <- c(OR11 - OR10, OR11 - OR01, OR11)
  sub_vcov <- vcov_mat[c(idx1, idx2, idx_int), c(idx1, idx2, idx_int)]
  var_RERI <- t(G) %*% sub_vcov %*% G
  se_RERI <- sqrt(var_RERI)

  RERI_lower <- RERI - 1.96 * se_RERI
  RERI_upper <- RERI + 1.96 * se_RERI
  z <- RERI / se_RERI
  p_RERI <- 2 * pnorm(-abs(z))

  p1_orig <- name_map[p1_valid]
  p2_orig <- name_map[p2_valid]

  return(data.frame(
    Pollutant_A = p1_orig,
    Pollutant_B = p2_orig,
    RERI = RERI,
    RERI_lower = RERI_lower,
    RERI_upper = RERI_upper,
    RERI_p = p_RERI,
    AP = AP,
    S = S,
    stringsAsFactors = FALSE
  ))
}

n <- length(pollutants_valid)
pair_indices <- combn(n, 2, simplify = FALSE)

results_list <- pblapply(pair_indices, function(idx) {
  p1 <- pollutants_valid[idx[1]]
  p2 <- pollutants_valid[idx[2]]
  calc_interaction(p1, p2, data_bin, "Hypertension", covariates_valid)
})

results_list <- results_list[!sapply(results_list, is.null)]

if (length(results_list) == 0) {
  stop("No valid pollutant pairs could be computed; check data and model convergence.")
}

results_df <- do.call(rbind, results_list)

# -----------------------------------------------------------------------------
# FDR correction and output
# -----------------------------------------------------------------------------
results_df$FDR_q <- p.adjust(results_df$RERI_p, method = "BH")
results_df$Significant <- ifelse(results_df$FDR_q < 0.05, "Yes", "No")

output <- results_df %>%
  mutate(
    RERI_CI = paste0(round(RERI, 3), " (", round(RERI_lower, 3), ", ",
                     round(RERI_upper, 3), ")"),
    RERI_Raw_P = format(RERI_p, scientific = FALSE, digits = 4),
    AP_Estimate = round(AP, 3),
    S_Estimate = round(S, 3),
    FDR_q = format(FDR_q, scientific = FALSE, digits = 4)
  ) %>%
  select(Pollutant_A, Pollutant_B, RERI_CI, RERI_Raw_P, AP_Estimate, S_Estimate,
         FDR_q, Significant)

output <- output[order(output$FDR_q), ]

cat("Computed", nrow(output), "pollutant pairs.\n")
print(head(output, 10))
