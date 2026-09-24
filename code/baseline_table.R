# =============================================================================
# Table S6: background information of study participants by hypertension status
# Three-line table, Times New Roman 12 pt
# Variables match the ML pipeline: covariates_cont + covariates_cat +
# seven fasting clinical indicators (PLT, WBC, Hb, ALB, GLO, UA, HCY)
# =============================================================================

rm(list = ls())

# Load packages
pacman::p_load(
  tidyverse, readxl, gtsummary, flextable, officer
)

# -----------------------------------------------------------------------------
# Input (load before running)
# -----------------------------------------------------------------------------
# df <- readxl::read_xlsx(...) %>% as_tibble()

# Outcome
df <- df %>%
  mutate(Hypertension = factor(Hypertension, levels = c(0, 1),
                               labels = c("Non-Hypertension", "Hypertension")))

# Categorical recoding
df <- df %>%
  mutate(
    Gender = factor(Gender, levels = c(1, 2),
                    labels = c("Male", "Female")),
    Smoke  = factor(Smoke, levels = c(1, 2, 3),
                    labels = c("Never", "Former", "Current")),
    Drink  = factor(Drink, levels = c(1, 2, 3, 4),
                    labels = c("Never", "Once in a while", "Often", "Every day")),
    Exercisefrequency = factor(Exercisefrequency, levels = c(1, 2, 3, 4),
                               labels = c("Every day", "Once a week",
                                          "Once in a while", "Never")),
    Eatinghabits = factor(Eatinghabits, levels = c(1, 2, 3),
                          labels = c("Equilibrium", "Vegetarian", "Non-vegetarian"))
  )

# -----------------------------------------------------------------------------
# Variable list
# -----------------------------------------------------------------------------
var_list <- list(
  # Baseline covariates (continuous)
  "Age2"                         = "Age (years)",
  "BMI"                          = "BMI (kg/m²)",
  "Leftfastingbloodglucose2mmol" = "Fasting Glucose (mmol/L)",
  # Blood pressure (not used in ML models)
  "LeftSBP2mmHg"                 = "SBP (mmHg)",
  "LeftDBP2mmHg"                 = "DBP (mmHg)",
  # Baseline covariates (categorical)
  "Gender"                       = "Gender",
  "Exercisefrequency"            = "Physical Activity",
  "Eatinghabits"                 = "Dietary Habit",
  "Smoke"                        = "Smoking Status",
  "Drink"                        = "Alcohol Drinking",
  # Seven fasting clinical indicators (ML predictors)
  "Platelet2gl"                  = "Platelet (10⁹/L)",
  "WBC2gl"                       = "WBC (10⁹/L)",
  "Hemoglobin2gL"                = "Hemoglobin (g/L)",
  "ALB"                          = "Albumin (g/L)",
  "GLO"                          = "Globulin (g/L)",
  "UA"                           = "Uric Acid (μmol/L)",
  "HCY"                          = "Homocysteine (μmol/L)"
)

# Keep only variables present in the data
available_vars <- intersect(names(var_list), colnames(df))
var_list <- var_list[available_vars]
cat("Variables included in the baseline table:", length(var_list), "\n")

# -----------------------------------------------------------------------------
# Analysis subset
# -----------------------------------------------------------------------------
df_analysis <- df %>% select(Hypertension, all_of(available_vars))

# -----------------------------------------------------------------------------
# gtsummary baseline table
# -----------------------------------------------------------------------------
N_total <- nrow(df_analysis)
n_non   <- sum(df_analysis$Hypertension == "Non-Hypertension")
n_ht    <- sum(df_analysis$Hypertension == "Hypertension")

tbl <- df_analysis %>%
  tbl_summary(
    by        = Hypertension,
    label     = var_list,
    statistic = list(
      all_continuous()  ~ "{mean} ± {sd}",
      all_categorical()  ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous()  ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    missing = "no"
  ) %>%
  add_overall(
    last      = FALSE,
    statistic = list(
      all_continuous()  ~ "{mean} ± {sd}",
      all_categorical()  ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous()  ~ 2,
      all_categorical() ~ c(0, 1)
    )
  ) %>%
  add_p(
    test = list(
      all_continuous()  ~ "t.test",
      all_categorical() ~ "fisher.test"
    ),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  modify_header(
    label   = "**Characteristic**",
    stat_0  = paste0("**Total**  \n(N = ", N_total, ")"),
    stat_1  = paste0("**Non-Hypertension**  \nN = ", n_non),
    stat_2  = paste0("**Hypertension**  \nN = ", n_ht),
    p.value = "**p-value**"
  ) %>%
  modify_footnote(everything() ~ NA) %>%
  bold_labels()

# -----------------------------------------------------------------------------
# Flextable (three-line table, Times New Roman 12 pt)
# -----------------------------------------------------------------------------
ft <- tbl %>%
  gtsummary::as_flex_table() %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::fontsize(size = 12, part = "all") %>%
  flextable::align(j = 1, align = "left",  part = "all") %>%
  flextable::align(j = 2:5, align = "center", part = "all") %>%
  flextable::align(align = "center", part = "header") %>%
  flextable::hline_top(border = officer::fp_border(width = 1.5), part = "header") %>%
  flextable::hline_bottom(border = officer::fp_border(width = 0.5), part = "header") %>%
  flextable::hline_bottom(border = officer::fp_border(width = 1.5), part = "body") %>%
  flextable::autofit() %>%
  flextable::padding(padding.top = 2, padding.bottom = 2, part = "all")

ft
