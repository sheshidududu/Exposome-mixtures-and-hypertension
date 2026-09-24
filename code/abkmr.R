# =============================================================================
# Adaptive Bayesian Kernel Machine Regression (aBKMR) pipeline
# Outcome: hypertension
# Notes:
#   - a_kmbayes() runs all "iter" iterations; bkmr post-processing functions
#     (ExtractPIPs, OverallRiskSummaries, etc.) automatically use the second
#     half of the iterations as post-burnin samples (burnin = iter / 2).
# =============================================================================

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 1. Load packages
# -----------------------------------------------------------------------------
pkg <- c("dplyr", "ggplot2", "aBKMR", "bkmr")
for (i in pkg) {
  if (!require(i, character.only = TRUE)) {
    stop(paste0("Package not installed: ", i))
  }
}

# -----------------------------------------------------------------------------
# 2. Global parameters
# -----------------------------------------------------------------------------
set.seed(123)

outcome_var <- "Hypertension"
iter <- 100000    # total MCMC iterations -> post-burnin = 50000
family <- "binomial"

# -----------------------------------------------------------------------------
# 3. Input (load before running)
# -----------------------------------------------------------------------------
# data: exposure data with the pollutants listed below, baseline covariates
# and the binary outcome "Hypertension"
# data <- read_excel(...)
data <- as.data.frame(data)

# -----------------------------------------------------------------------------
# 4. Create combined exposures for highly collinear pairs
# -----------------------------------------------------------------------------
data$DCCA <- rowMeans(data[, c("trans-DCCA", "cis-DCCA")], na.rm = TRUE)
data$TETDOXY <- rowMeans(data[, c("TET", "DOXY")], na.rm = TRUE)

# -----------------------------------------------------------------------------
# 5. Exposure list (66 individual pollutants + 2 combined = 68)
# -----------------------------------------------------------------------------
pollutants <- c(
  "SFM", "4AN6C2B", "SFSX", "TMTH", "SFTHI",
  "NOR", "ENOR", "CIP", "OFLX", "OTC", "AM",
  "IMI", "NIT", "ACE", "THM", "CLO", "THI",
  "IPP", "CYC", "DIN", "IMID", "DM-THM",
  "DN-IMI", "5-OH-IMI", "DM-ACE",
  "FP", "DMM", "CBDZ", "CPS-O",
  "DZN", "DDVP", "MALT", "3OH-CBF",
  "ATZ", "ATZ-DE-DSI-2HD",
  "BPF", "BPS", "BPC", "BPE", "BPP",
  "2-OHN", "1-OHN",
  "3-OHF", "2-OHF",
  "2-OHP", "1-OHP",
  "4-OHP", "1-OHPYR",
  "3-PBA", "4-F-3-PBA",
  "2-MB-3CA",
  "IMI-olefin",
  "3,5,6-TCP",
  "V", "Fe", "Co", "Cu", "Zn",
  "As", "Se", "Sr", "Mo", "Cd",
  "Te", "Tl", "Pb",
  "DCCA",
  "TETDOXY"
)

# -----------------------------------------------------------------------------
# 6. Covariates
# -----------------------------------------------------------------------------
cov_cat <- c("Gender", "Exercisefrequency", "Eatinghabits", "Smoke", "Drink")
cov_con <- c("Leftfastingbloodglucose2mmol", "BMI", "Age2")

# -----------------------------------------------------------------------------
# 7. Check variables
# -----------------------------------------------------------------------------
need_var <- c(outcome_var, pollutants, cov_cat, cov_con)
miss_var <- need_var[!need_var %in% names(data)]
if (length(miss_var) > 0) {
  stop(paste("Missing variables:", paste(miss_var, collapse = ", ")))
}

# -----------------------------------------------------------------------------
# 8. Convert factors
# -----------------------------------------------------------------------------
data[cov_cat] <- lapply(data[cov_cat], factor)

# -----------------------------------------------------------------------------
# 9. Convert outcome to 0/1
# -----------------------------------------------------------------------------
if (is.factor(data[[outcome_var]])) {
  if (nlevels(data[[outcome_var]]) != 2) stop("Outcome must be binary.")
  y <- as.numeric(data[[outcome_var]]) - 1
} else {
  y <- as.numeric(data[[outcome_var]])
}
if (!all(y %in% c(0, 1))) stop("Outcome is not coded as 0/1.")

# -----------------------------------------------------------------------------
# 10. Missing check
# -----------------------------------------------------------------------------
analysis_data <- data[, c(outcome_var, pollutants, cov_cat, cov_con)]

na_summary <- data.frame(
  Variable = names(analysis_data),
  Missing = sapply(analysis_data, function(x) sum(is.na(x)))
)
na_summary$Percent <- round(100 * na_summary$Missing / nrow(analysis_data), 2)

# -----------------------------------------------------------------------------
# 11. Complete cases
# -----------------------------------------------------------------------------
complete_index <- complete.cases(analysis_data)
cat("Original sample :", nrow(data), "\n")
cat("Complete sample :", sum(complete_index), "\n")

data <- data[complete_index, ]
y <- y[complete_index]

# -----------------------------------------------------------------------------
# 12. Variable mapping (original <-> R-safe names)
# -----------------------------------------------------------------------------
map <- data.frame(
  Original = pollutants,
  R_Name = make.names(pollutants, unique = TRUE)
)

# -----------------------------------------------------------------------------
# 13. Create exposure matrix
# -----------------------------------------------------------------------------
Z.raw <- as.matrix(data[, pollutants])
colnames(Z.raw) <- make.names(colnames(Z.raw), unique = TRUE)

Z <- scale(Z.raw)
storage.mode(Z) <- "double"

# -----------------------------------------------------------------------------
# 14. Covariate matrix
# -----------------------------------------------------------------------------
X <- model.matrix(
  ~ Gender + Exercisefrequency + Eatinghabits + Smoke + Drink +
    Leftfastingbloodglucose2mmol + BMI + Age2,
  data = data
)[, -1]
storage.mode(X) <- "double"

# -----------------------------------------------------------------------------
# 15. Final data check
# -----------------------------------------------------------------------------
stopifnot(nrow(Z) == length(y))
stopifnot(nrow(X) == length(y))
stopifnot(nrow(Z) == nrow(X))

storage.mode(y) <- "double"

# -----------------------------------------------------------------------------
# 16. Select representative knots
# -----------------------------------------------------------------------------
k_id <- sam_py_r(R = Z, nd = 100, num_nn = 100, w = FALSE)

# -----------------------------------------------------------------------------
# 17. Model fitting
# -----------------------------------------------------------------------------
start_time <- Sys.time()

fit_abkmr <- a_kmbayes(
  y = y,
  Z = Z,
  X = X,
  knots = Z[k_id, ],
  family = family,
  iter = iter,
  varsel = TRUE,
  est.h = TRUE,
  verbose = TRUE
)

end_time <- Sys.time()
runtime <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)

cat("Model finished.\n")
cat("Running time:", runtime, "minutes\n")

# -----------------------------------------------------------------------------
# 18. Posterior inclusion probability (PIP)
# -----------------------------------------------------------------------------
pip <- ExtractPIPs(fit_abkmr)

if (nrow(pip) != ncol(Z)) stop("PIP extraction failed.")

# Restore original variable names
pip <- pip |>
  left_join(map, by = c("variable" = "R_Name"))

pip$Display_Name <- ifelse(is.na(pip$Original), pip$variable, pip$Original)

pip <- pip |>
  arrange(desc(PIP))

rownames(pip) <- NULL
pip$Rank <- 1:nrow(pip)

# -----------------------------------------------------------------------------
# 19. PIP figure (top 30)
# -----------------------------------------------------------------------------
plot_pip <- pip[1:30, ]
plot_pip$Display_Name <- factor(plot_pip$Display_Name,
                                levels = rev(plot_pip$Display_Name))

p_pip <- ggplot(plot_pip, aes(Display_Name, PIP)) +
  geom_col(fill = "#2C7FB8", width = .75) +
  geom_text(aes(label = sprintf("%.3f", PIP)), hjust = -0.15, size = 3.5) +
  coord_flip() +
  expand_limits(y = 1.05) +
  labs(x = NULL, y = "Posterior Inclusion Probability") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    plot.title = element_text(face = "bold", hjust = .5)
  )

ggsave("Figure1_PIP.pdf", plot = p_pip, width = 8, height = 10)
ggsave("Figure1_PIP.png", plot = p_pip, dpi = 600, width = 8, height = 10)

# -----------------------------------------------------------------------------
# 20. Overall mixture risk
# -----------------------------------------------------------------------------
quants <- seq(0.25, 0.75, by = 0.05)

newz_list <- lapply(1:ncol(Z), function(i) {
  data.frame(x = quantile(Z[, i], probs = quants, na.rm = TRUE))
})
newz <- do.call(cbind, newz_list)
colnames(newz) <- colnames(Z)

overall_res <- a_OverallRiskSummaries_vary(
  fit_abkmr,
  newz = newz,
  data.comps = fit_abkmr$data.comps
)

risk_df <- overall_res$risk_overall
risk_df$quantile <- quants

p_overall <- ggplot(risk_df, aes(x = quantile, y = est_p)) +
  geom_errorbar(aes(ymin = est_p - 1.96 * sd_p, ymax = est_p + 1.96 * sd_p),
                width = 0.02, color = "black") +
  geom_point(size = 2, color = "black") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  labs(x = "Mixture Quantile", y = "Overall Log-Odds Ratio") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"))

ggsave("Figure2_OverallRisk.pdf", p_overall, width = 7, height = 5)
ggsave("Figure2_OverallRisk.png", p_overall, width = 7, height = 5, dpi = 600)

# -----------------------------------------------------------------------------
# 21. Univariate exposure-response
# -----------------------------------------------------------------------------
single_res <- a_PredictorResponseUnivar(
  fit = fit_abkmr,
  quants = seq(0.1, 0.9, by = 0.1),
  method = "approx",
  data.comps = fit_abkmr$data.comps
)

single_plot <- dplyr::bind_rows(
  lapply(single_res, function(x) x$risk_overall)
)

if (nrow(single_plot) == 0) stop("No univariate results extracted.")

# Top 10 pollutants by PIP
top_vars <- pip$variable[1:10]

single_plot <- single_plot |>
  dplyr::filter(vars %in% top_vars)

# Restore original names
single_plot <- single_plot |>
  dplyr::left_join(map, by = c("vars" = "R_Name"))

single_plot$Display_Name <- ifelse(is.na(single_plot$Original),
                                   single_plot$vars, single_plot$Original)

p_single <- ggplot(single_plot, aes(x = z, y = est)) +
  geom_ribbon(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se), alpha = .20) +
  geom_line(linewidth = .9) +
  facet_wrap(~ Display_Name, scales = "free_x", ncol = 5) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  ) +
  labs(x = "Standardized Exposure", y = "Estimated Effect")

ggsave("Figure3_Univariate.pdf", p_single, width = 14, height = 8)
ggsave("Figure3_Univariate.png", p_single, width = 14, height = 8, dpi = 600)

# -----------------------------------------------------------------------------
# 22. Interaction analysis
# -----------------------------------------------------------------------------
inter_res <- a_PredictorResponseBivar(
  fit = fit_abkmr,
  data.comps = fit_abkmr$data.comps,
  min.plot.dist = 1
)

inter_plot <- PredictorResponseBivarLevels(
  inter_res,
  Z = fit_abkmr$Z,
  qs = c(0.25, 0.50, 0.75)
)

# Keep top 5 pollutants by PIP
top5 <- pip$variable[1:5]

inter_plot <- inter_plot |>
  dplyr::filter(variable1 %in% top5, variable2 %in% top5) |>
  dplyr::filter(variable1 != variable2)

# Recover original names
inter_plot <- inter_plot |>
  left_join(map, by = c("variable1" = "R_Name")) |>
  rename(variable1_name = Original)

inter_plot <- inter_plot |>
  left_join(map, by = c("variable2" = "R_Name")) |>
  rename(variable2_name = Original)

inter_plot$variable1_name <- ifelse(is.na(inter_plot$variable1_name),
                                    inter_plot$variable1, inter_plot$variable1_name)
inter_plot$variable2_name <- ifelse(is.na(inter_plot$variable2_name),
                                    inter_plot$variable2, inter_plot$variable2_name)

inter_plot_clean <- inter_plot %>% filter(!is.na(est))

p_inter <- ggplot(inter_plot_clean, aes(x = z1, y = est, group = quantile, colour = quantile)) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  facet_grid(variable2_name ~ variable1_name, scales = "free_x") +
  scale_colour_manual(values = c("0.25" = "#1B9E77", "0.5" = "#D95F02", "0.75" = "#7570B3")) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black")) +
  labs(x = "Standardized exposure", y = "Estimated effect", colour = "Modifier quantile")

ggsave("Figure4_Interaction.pdf", p_inter, width = 12, height = 10)
ggsave("Figure4_Interaction.png", p_inter, width = 12, height = 10, dpi = 600)

# -----------------------------------------------------------------------------
# 23. Three-level forest plot (top 10 and top 20)
# -----------------------------------------------------------------------------
point25 <- apply(Z, 2, quantile, probs = 0.25, na.rm = TRUE)
point50 <- apply(Z, 2, quantile, probs = 0.50, na.rm = TRUE)
point75 <- apply(Z, 2, quantile, probs = 0.75, na.rm = TRUE)

res25 <- a_PredictorResponseUnivar(
  fit = fit_abkmr,
  quants = seq(0.10, 0.90, 0.10),
  point1 = point25,
  data.comps = fit_abkmr$data.comps
)
res50 <- a_PredictorResponseUnivar(
  fit = fit_abkmr,
  quants = seq(0.10, 0.90, 0.10),
  point1 = point50,
  data.comps = fit_abkmr$data.comps
)
res75 <- a_PredictorResponseUnivar(
  fit = fit_abkmr,
  quants = seq(0.10, 0.90, 0.10),
  point1 = point75,
  data.comps = fit_abkmr$data.comps
)

extract_est <- function(res, level) {
  do.call(rbind, lapply(res, function(x) {
    data.frame(
      variable = x$risk_overall$vars[1],
      Estimate = x$est$est,
      SE = x$est$sd,
      Level = level,
      stringsAsFactors = FALSE
    )
  }))
}

forest25 <- extract_est(res25, "25%")
forest50 <- extract_est(res50, "50%")
forest75 <- extract_est(res75, "75%")

forest <- rbind(forest25, forest50, forest75)

forest$Lower <- forest$Estimate - 1.96 * forest$SE
forest$Upper <- forest$Estimate + 1.96 * forest$SE

# Restore variable names
forest$variable <- map$Original[match(forest$variable, map$R_Name)]

draw_forest <- function(forest, top_vars, title, filename, height) {
  forest_plot <- subset(forest, variable %in% top_vars)
  forest_plot$variable <- factor(forest_plot$variable, levels = rev(top_vars))

  forest_plot$y_base <- as.numeric(forest_plot$variable)
  forest_plot$y <- forest_plot$y_base +
    ifelse(forest_plot$Level == "25%", -0.25,
           ifelse(forest_plot$Level == "50%", 0, 0.25))

  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60", linewidth = 0.5) +
    geom_segment(data = forest_plot,
                 aes(x = Lower, xend = Upper, y = y, yend = y, colour = Level),
                 linewidth = 0.8) +
    geom_segment(data = forest_plot,
                 aes(x = Estimate, xend = Estimate,
                     y = y - 0.08, yend = y + 0.08, colour = Level),
                 linewidth = 0.8) +
    geom_point(data = forest_plot,
               aes(x = Estimate, y = y, colour = Level),
               shape = 18, size = 4.5) +
    scale_y_continuous(breaks = unique(forest_plot$y_base), labels = rev(top_vars)) +
    scale_colour_manual(
      values = c("25%" = "#2ECC71", "50%" = "#1A5276", "75%" = "#C0392B"),
      name = "Overall mixture"
    ) +
    labs(x = "Estimated effect (95% CI)", y = NULL, title = title) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "right",
      axis.text.y = element_text(face = "bold", size = 11),
      axis.text.x = element_text(size = 11),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )

  ggsave(paste0(filename, ".pdf"), p, width = 8, height = height)
  ggsave(paste0(filename, ".png"), p, width = 8, height = height, dpi = 600)
}

top10_vars <- map$Original[match(pip$variable[1:10], map$R_Name)]
draw_forest(forest, top10_vars, "Top 10 Variables", "Figure5_Forest_Top10", 6)

top20_vars <- map$Original[match(pip$variable[1:20], map$R_Name)]
draw_forest(forest, top20_vars, "Top 20 Variables", "Figure5_Forest_Top20", 10)

# -----------------------------------------------------------------------------
# 24. MCMC convergence diagnostics (post-burnin samples only)
# -----------------------------------------------------------------------------
library(coda)

n_iter <- nrow(fit_abkmr$r)
burnin_end <- floor(n_iter / 2)
post_idx <- (burnin_end + 1):n_iter

cat("Total iterations:", n_iter, "\n")
cat("Burn-in discarded:", burnin_end, "\n")
cat("Post-burnin samples:", length(post_idx), "\n")

# Select top-3 exposures by PIP and first 2 covariates
exposure_names <- colnames(fit_abkmr$Z)
n_exp <- length(exposure_names)

cov_names <- colnames(fit_abkmr$X)
n_cov <- length(cov_names)

if (exists("pip") && nrow(pip) > 0) {
  pip_sorted <- pip %>% arrange(desc(PIP))
  top_exp_rnames <- pip_sorted$variable[1:min(3, nrow(pip_sorted))]
  top_exp_idx <- match(top_exp_rnames, exposure_names)
  top_exp_idx <- top_exp_idx[!is.na(top_exp_idx)]
  top_exp_rnames <- exposure_names[top_exp_idx]

  if (length(top_exp_idx) < 3) {
    add_idx <- setdiff(1:3, top_exp_idx)
    add_idx <- add_idx[add_idx <= n_exp]
    top_exp_idx <- unique(c(top_exp_idx, add_idx))
    top_exp_rnames <- exposure_names[top_exp_idx]
  }
} else {
  top_exp_idx <- 1:min(3, n_exp)
  top_exp_rnames <- exposure_names[top_exp_idx]
}

cov_idx <- if (n_cov >= 2) 1:2 else 1:n_cov
selected_cov_names <- cov_names[cov_idx]

# Acceptance rate
acc_result <- data.frame(
  Parameter = c("r (average)", "lambda"),
  Mean = c(
    mean(fit_abkmr$acc.r[post_idx], na.rm = TRUE),
    mean(fit_abkmr$acc.lambda[post_idx], na.rm = TRUE)
  )
)
cat("Acceptance rate (post-burnin):\n")
print(acc_result)

# Trace plots (2x2 layout)
par(mfrow = c(2, 2))

plot(post_idx, as.numeric(fit_abkmr$lambda)[post_idx], type = "l",
     main = "Lambda (Variable Selection Penalty) - Post-Burnin",
     xlab = "Iteration", ylab = "Value", col = "#1f78b4")

if (length(top_exp_idx) >= 1) {
  plot(post_idx, fit_abkmr$r[post_idx, top_exp_idx[1]], type = "l",
       main = paste("r (Inclusion) for", top_exp_rnames[1]),
       xlab = "Iteration", ylab = "Value", col = "#33a02c")
} else {
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(1, 1, "No exposure selected")
}

if (length(top_exp_idx) >= 2) {
  plot(post_idx, fit_abkmr$r[post_idx, top_exp_idx[2]], type = "l",
       main = paste("r (Inclusion) for", top_exp_rnames[2]),
       xlab = "Iteration", ylab = "Value", col = "#33a02c")
} else {
  plot(post_idx, as.numeric(fit_abkmr$lambda)[post_idx], type = "l",
       main = "Lambda (Fallback)", xlab = "Iteration", ylab = "Value", col = "#1f78b4")
}

if (length(cov_idx) >= 1) {
  plot(post_idx, fit_abkmr$beta[post_idx, cov_idx[1]], type = "l",
       main = paste("Beta for", selected_cov_names[1]),
       xlab = "Iteration", ylab = "Value", col = "#e31a1c")
} else {
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(1, 1, "No covariate selected")
}

# Posterior density (2x2 layout)
par(mfrow = c(2, 2))

plot(density(as.numeric(fit_abkmr$lambda)[post_idx], na.rm = TRUE),
     main = "Lambda - Post-Burnin", col = "#1f78b4", lwd = 2)

if (length(top_exp_idx) >= 1) {
  plot(density(fit_abkmr$r[post_idx, top_exp_idx[1]], na.rm = TRUE),
       main = paste("r for", top_exp_rnames[1]), col = "#33a02c", lwd = 2)
} else {
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(1, 1, "No exposure selected")
}

if (length(top_exp_idx) >= 2) {
  plot(density(fit_abkmr$r[post_idx, top_exp_idx[2]], na.rm = TRUE),
       main = paste("r for", top_exp_rnames[2]), col = "#33a02c", lwd = 2)
} else {
  plot(density(as.numeric(fit_abkmr$lambda)[post_idx], na.rm = TRUE),
       main = "Lambda (Fallback)", col = "#1f78b4", lwd = 2)
}

if (length(cov_idx) >= 1) {
  plot(density(fit_abkmr$beta[post_idx, cov_idx[1]], na.rm = TRUE),
       main = paste("Beta for", selected_cov_names[1]), col = "#e31a1c", lwd = 2)
} else {
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(1, 1, "No covariate selected")
}

# Autocorrelation (2x2 layout)
par(mfrow = c(2, 2))

acf(as.numeric(fit_abkmr$lambda)[post_idx],
    main = "Lambda - Post-Burnin", col = "#1f78b4", lwd = 2)

if (length(top_exp_idx) >= 1) {
  acf(fit_abkmr$r[post_idx, top_exp_idx[1]],
      main = paste("r for", top_exp_rnames[1]), col = "#33a02c", lwd = 2)
} else {
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(1, 1, "No exposure selected")
}

if (length(top_exp_idx) >= 2) {
  acf(fit_abkmr$r[post_idx, top_exp_idx[2]],
      main = paste("r for", top_exp_rnames[2]), col = "#33a02c", lwd = 2)
} else {
  acf(as.numeric(fit_abkmr$lambda)[post_idx],
      main = "Lambda (Fallback)", col = "#1f78b4", lwd = 2)
}

if (length(cov_idx) >= 1) {
  acf(fit_abkmr$beta[post_idx, cov_idx[1]],
      main = paste("Beta for", selected_cov_names[1]), col = "#e31a1c", lwd = 2)
} else {
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(1, 1, "No covariate selected")
}

# Effective sample size
cat("\nComputing effective sample size (ESS, post-burnin only)...\n")

ess_list <- list()

ess_list[["Lambda"]] <- effectiveSize(
  as.mcmc(as.numeric(fit_abkmr$lambda)[post_idx])
)

for (i in seq_along(top_exp_idx)) {
  param_name <- paste0("r_", top_exp_rnames[i])
  ess_list[[param_name]] <- effectiveSize(
    as.mcmc(fit_abkmr$r[post_idx, top_exp_idx[i]])
  )
}

for (i in seq_along(cov_idx)) {
  param_name <- paste0("beta_", selected_cov_names[i])
  ess_list[[param_name]] <- effectiveSize(
    as.mcmc(fit_abkmr$beta[post_idx, cov_idx[i]])
  )
}

ess_df <- data.frame(
  Parameter = names(ess_list),
  ESS = round(unlist(ess_list), 2)
)

ess_df <- ess_df %>% arrange(ESS)

cat("ESS summary (post-burnin key parameters):\n")
print(ess_df)

cat("\naBKMR analysis pipeline complete.\n")
