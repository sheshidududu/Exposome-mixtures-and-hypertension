# =============================================================================
# Mann-Whitney U tests: pollutant concentrations by hypertension status,
# with per-group medians and 95% CIs (percentile bootstrap)
# =============================================================================

library(readxl)
library(dplyr)
library(purrr)
library(writexl)
library(boot)

# -----------------------------------------------------------------------------
# Bootstrap median confidence interval
# -----------------------------------------------------------------------------
median_ci <- function(x, conf = 0.95, R = 1000) {
  if (length(x) < 3) {
    return(c(median(x, na.rm = TRUE), NA, NA))
  }
  boot_median <- function(data, indices) {
    median(data[indices], na.rm = TRUE)
  }
  boot_obj <- boot(x, boot_median, R = R)
  ci <- boot.ci(boot_obj, type = "perc", conf = conf)
  lower <- ci$perc[4]
  upper <- ci$perc[5]
  return(c(median(x, na.rm = TRUE), lower, upper))
}

# -----------------------------------------------------------------------------
# Input (load before running)
# -----------------------------------------------------------------------------
# data <- read_excel(...)

# -----------------------------------------------------------------------------
# Pollutants (70)
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

# -----------------------------------------------------------------------------
# Preprocessing
# -----------------------------------------------------------------------------
data$Hypertension <- factor(data$Hypertension, levels = c(0, 1),
                            labels = c("Non-HTN", "HTN"))
data <- data %>% mutate(across(all_of(pollutants), as.numeric))

# -----------------------------------------------------------------------------
# Mann-Whitney U tests + medians with bootstrap CIs
# -----------------------------------------------------------------------------
results <- map_df(pollutants, function(var) {
  group0 <- data[[var]][data$Hypertension == "Non-HTN"]
  group1 <- data[[var]][data$Hypertension == "HTN"]
  group_all <- data[[var]][!is.na(data[[var]])]

  group0 <- group0[!is.na(group0)]
  group1 <- group1[!is.na(group1)]

  n1 <- length(group0)
  n2 <- length(group1)

  ci_all <- median_ci(group_all, R = 1000)
  med_all <- ci_all[1]
  low_all <- ci_all[2]
  high_all <- ci_all[3]

  total_ci_str <- ifelse(!is.na(low_all) & !is.na(high_all),
                         sprintf("%.3f (%.3f-%.3f)", med_all, low_all, high_all),
                         sprintf("%.3f (NA-NA)", med_all))

  if (n1 > 0 & n2 > 0) {
    test <- wilcox.test(group0, group1, exact = FALSE, conf.int = FALSE)
    U <- test$statistic
    p <- test$p.value

    mu <- n1 * n2 / 2
    sigma <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    Z <- (U - mu) / sigma

    ci0 <- median_ci(group0, R = 1000)
    med0 <- ci0[1]
    low0 <- ci0[2]
    high0 <- ci0[3]

    nonhtn_ci_str <- ifelse(!is.na(low0) & !is.na(high0),
                            sprintf("%.3f (%.3f-%.3f)", med0, low0, high0),
                            sprintf("%.3f (NA-NA)", med0))

    ci1 <- median_ci(group1, R = 1000)
    med1 <- ci1[1]
    low1 <- ci1[2]
    high1 <- ci1[3]

    htn_ci_str <- ifelse(!is.na(low1) & !is.na(high1),
                         sprintf("%.3f (%.3f-%.3f)", med1, low1, high1),
                         sprintf("%.3f (NA-NA)", med1))

    tibble(
      Pollutant = var,
      `U statistic` = as.numeric(U),
      `Z-value` = Z,
      `P-value` = p,
      `Median (95% CI) - Non-HTN` = nonhtn_ci_str,
      `Median (95% CI) - HTN` = htn_ci_str,
      `Median (95% CI) - Total` = total_ci_str
    )
  } else {
    tibble(
      Pollutant = var,
      `U statistic` = NA_real_,
      `Z-value` = NA_real_,
      `P-value` = NA_real_,
      `Median (95% CI) - Non-HTN` = NA_character_,
      `Median (95% CI) - HTN` = NA_character_,
      `Median (95% CI) - Total` = total_ci_str
    )
  }
})

# -----------------------------------------------------------------------------
# FDR correction and rounding
# -----------------------------------------------------------------------------
results$`Adjusted P-value (FDR)` <- p.adjust(results$`P-value`, method = "fdr")

results <- results %>%
  mutate(
    `U statistic` = round(`U statistic`, 3),
    `Z-value` = round(`Z-value`, 3),
    `P-value` = round(`P-value`, 3),
    `Adjusted P-value (FDR)` = round(`Adjusted P-value (FDR)`, 3)
  )

results <- results %>%
  select(
    Pollutant,
    `U statistic`,
    `Z-value`,
    `P-value`,
    `Adjusted P-value (FDR)`,
    `Median (95% CI) - Non-HTN`,
    `Median (95% CI) - HTN`,
    `Median (95% CI) - Total`
  )

print(head(results, 10))
