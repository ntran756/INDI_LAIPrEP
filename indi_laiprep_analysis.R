# Reproducible code for: Intersectional discrimination and preferences for long-acting injectable HIV PrEP 
# among Black and Latino gay, bisexual, and men who have sex with men in Philadelphia: a cross-sectional study

# Set up ----
## load packages
library(tidyverse)
library(here)
library(marginaleffects)
library(tableone)

## import data
df <- read_csv(here("final_dta_20220327.csv")) |> 
  janitor::clean_names() |> 
  # filter to only cisgender men
  filter(gender == "Man") 

# Recode socio-demographics for analysis ----
df <- df |> 
  mutate(
    race_cat = factor(case_when(
      race == "Multiracial" & hisp_yn == "Yes" ~ 5, 
      race == "White" & hisp_yn == "No" ~ 4,
      race == "White" & hisp_yn == "Yes" ~ 3,
      str_detect(race, "Black") & hisp_yn == "No" ~ 2,
      str_detect(race, "Black") & hisp_yn == "Yes" ~ 1
    ), labels = c("Black Latino", "Black non-Latino", "White Latino", "White non-Latino", "multiracial Latino")),
    age_grp = factor(case_when(
      age < 25 ~ 1, 
      age >= 25 & age < 29 ~ 2, 
      age >= 30 & age < 34 ~ 3,
      T ~ 4
    ), labels = c("18-24", "25-29", "30-34", "35+")), 
    educ_3cat = factor(case_when(
      str_detect(educ, "High school|high school") ~ 1, 
      str_detect(educ, "2-year degree|Some college") ~ 2, 
      T ~ 3
    ), labels = c("High school or less", "Some college", "4-year degree or higher")),
    employed = factor(case_when(
      is.na(employ) ~ NA,
      str_detect(employ, "No") ~ 1, 
      str_detect(employ, "Part-time") ~ 2, 
      T ~ 3
    ), labels = c("Unemployed", "Part-time", "Full-time")),
    income_3cat = factor(case_when(
      income %in% c("$0 - $14,999", "$15,000 - $24,999", "$25,000 - $34,999") ~ 1, 
      income %in% c("$35,000 - $44,999", "$45,000 - $54,999") ~ 2,
      T ~ 3
    ), labels = c("$0-$34,999", "$35,000-$54,999", "$55,000+")),
    unstable_house = factor(
      ifelse(str_detect(housing, "Homeless|Temporarily"), 1, 0), 
      labels = c("No", "Yes")
    ),
    no_insure = factor(
      ifelse(str_detect(insure, "I do not have"), 1, 0), 
      labels = c("No", "Yes")
    )
  )

# Calculate intersectional anticipated, daily enacted, and major enacted discrimination
df <- df |> 
  mutate(
    across(
      in_di_a_1:in_di_a_9, 
      ~ case_when(
        .x == "Strongly agree" ~ 4, 
        .x == "Agree" ~ 3,
        .x == "Neither agree nor disagree" ~ 2,
        .x == "Disagree" ~ 1, 
        .x == "Strongly disagree" ~ 0,
      )
    ), 
    across(
      in_di_d_1:in_di_d_9, 
      ~ case_when(
        .x %in% c("Never", "Yes, but not in the past year") ~ 0,
        .x == "Yes, once or twice in the past year" ~ 1,
        .x == "Yes, many times in the past year" ~ 2
      )
    ),
    across(
      in_di_m_1:in_di_m_13, 
      ~ case_when(
        .x == "Never" ~ 0, 
        .x == "Once" ~ 1, 
        .x == "More than Once" ~ 2
      )
    ), 
    # calculate scores
    indi_a = rowMeans(across(in_di_a_1:in_di_a_9), na.rm = T),
    indi_d = rowSums(across(in_di_d_1:in_di_d_9), na.rm = T),
    indi_m = rowSums(across(in_di_m_1:in_di_m_13), na.rm = T), 
    # group scores into percentiles
    indi_a_cat = cut(indi_a, breaks = quantile(indi_a, seq(0, 1, 0.25)), include.lowest = T),
    indi_d_cat = cut(indi_d, breaks = quantile(indi_d, seq(0, 1, 0.25)), include.lowest = T),
    indi_m_cat = cut(indi_m, breaks = quantile(indi_m, seq(0, 1, 0.25)), include.lowest = T),
    # recode percentiles as continuous measure to test linear trends
    indi_a_cont = as.numeric(indi_a_cat),
    indi_d_cont = as.numeric(indi_d_cat),
    indi_m_cont = as.numeric(indi_m_cat)
  ) 

# check distribution of discrimination scores
df |> group_by(indi_a_cat) |> summarise(n = n(), min(indi_a), max(indi_a))
df |> group_by(indi_d_cat) |> summarise(n = n(), min(indi_d), max(indi_d))
df |> group_by(indi_m_cat) |> summarise(n = n(), min(indi_m), max(indi_m))

# Evaluate reliability of each INDI measure
df_indi_a <- df |> select(in_di_a_1:in_di_a_9)
df_indi_d <- df |> select(in_di_d_1:in_di_d_9)
df_indi_m <- df |> select(in_di_m_1:in_di_m_13)

# calculate cronbach alpha 
ltm::cronbach.alpha(df_indi_a, B = 1000, CI = T, nr.rm = T)
ltm::cronbach.alpha(df_indi_d, B = 1000, CI = T, na.rm = T)
ltm::cronbach.alpha(df_indi_m, B = 1000, CI = T, na.rm = T)
rm(df_indi_a, df_indi_d, df_indi_m)

# Check prep use history
count(df, ever_prep, current_prep)

# any oral prep use: 1=yes, 0=no 
df <- df |> 
  mutate(no_ever_prep = ifelse(ever_prep == "No", 1, 0))

# prep knowledge
df <- df |> 
  mutate(
    # correct responses if selected "TRUE"
    across(
      c(prep_know_1:prep_know_4, prep_know_6, prep_know_8, prep_know_9),
      ~ ifelse(.x == "TRUE", 1, 0)
    ),
    # correct responses if selected "FALSE"
    across(
      c(prep_know_5, prep_know_7), 
      ~ ifelse(.x == "FALSE", 1, 0)
    ),
    prep_know_sum = rowSums(across(prep_know_1:prep_know_9), na.rm = T)
  )

# prefers injectable: 1=yes, 0=no
df <- df |> 
  mutate(
    prefer_inj_to_daily = ifelse(str_detect(injectable_willing_1, "Extremely likely|Somewhat likely"), 1, 0),
    prefer_inj_to_event = ifelse(str_detect(injectable_willing_2, "Extremely likely|Somewhat likely"), 1, 0),
    prefer_inj_to_daily_fct = factor(prefer_inj_to_daily, labels = c("No", "Yes")),
    prefer_inj_to_event_fct = factor(prefer_inj_to_event, labels = c("No", "Yes"))
  ) 

# check distribution 
count(df, prefer_inj_to_daily) |> mutate(pct = n/sum(n)*100)
count(df, prefer_inj_to_event) |> mutate(pct = n/sum(n)*100)

# Socio-demographic table 1 ----
v_names <- df |> 
  select(age_grp, age, sexor, educ_3cat, employed, income_3cat, 
         unstable_house, no_insure, prep_know_sum, ever_prep, 
         indi_a_cat, indi_a, indi_d_cat, indi_d, indi_m_cat, indi_m, 
         prefer_inj_to_daily, prefer_inj_to_event) |> 
  names()

v_fct_names <- discard(v_names, ~ .x %in% c("age", "indi_a", "indi_d", "indi_m"))

df_tab <- CreateTableOne(
  vars = v_names, 
  strata = "race_cat",
  data = df, 
  factorVars = v_fct_names, 
  includeNA = T, 
  addOverall = T
)

tab1 <- print(
  df_tab, 
  nonnormal = "age",
  quote = T,
  noSpaces = T, 
  printToggle = F
)

# Save table
# write.csv(tab1, here("output", "table1-participants_characteristics.csv"))
rm(df_tab, v_names, v_fct_names)

# calculate the number (%) in each indi quartile ----
tab2a <- df |> 
  pivot_longer(cols = c("indi_a_cat", "indi_d_cat", "indi_m_cat")) |>
  group_by(name, value) |> 
  summarise(
    n_total = n(),
    n_prefer = sum(prefer_inj_to_event == 1, na.rm = T),
    pct_prefer = round(n_prefer / n_total * 100, 1),
  )
tab2a

# save table
# write.csv(tab2a, here("output", "table2a-indi_quartiles.csv"))

# Regression ----
# helpful function to iterate the fitting multiple models for main analysis 
fit_model <- function(outcome, exposure, covariates = NULL, model_num, data) {
  predictors <- c(exposure, covariates)
  
  formula <- as.formula(
    paste(outcome, "~", paste(predictors, collapse = " + "))
  )
  
  model <- glm(formula, data = data, family = binomial())
  
  avg_comparisons(
    model,
    variables = exposure, 
    vcov = "HC3",
    comparison = "lnratioavg",
    transform = "exp"
  ) |>
    mutate(
      outcome = outcomes[[outcome]],
      model = paste0("m", model_num),
      exposure = outcomes[[exposure]],
      covariates = if (is.null(covariates)) "none" else paste(covariates, collapse = ", ")
    )
}

## specify variable for analysis
### exposures
exposure_var <- list(
  indi_a_cat = "InDI-A",
  indi_d_cat = "InDI-D",
  indi_m_cat = "InDI-M",
  indi_a_cont = "InDI-A_cont",
  indi_d_cont = "InDI-D_cont",
  indi_m_cont = "InDI-M_cont"
)

### define health outcomes and their display names
outcomes <- list(
  prefer_inj_to_daily = "Prefers Injectable to Daily Oral PrEP",
  prefer_inj_to_event = "Prefers Injectable to Event-Driven PrEP"
)
    
### define your models with different covariate combinations
model_specs <- list(
  list(covariates = NULL),
  list(covariates = c("age", "I(age^2)", "educ_3cat", "employed", "income_3cat", 
                      "sexor", "unstable_house", "race_cat", "prep_know_sum"))
)

### loop through functions to fit models and get outputs
df_reg <- bind_rows(
  lapply(names(outcomes), function(outcome_var) {
    lapply(names(exposure_var), function(expo_var) {
      lapply(seq_along(model_specs), function(model_num) {
        fit_model(
          outcome = outcome_var,
          exposure = expo_var,
          covariates = model_specs[[model_num]]$covariates,
          model_num = model_num,
          data = df)
      })
    })
  }) |> unlist(recursive = F)
)

### format outputs 
df_reg <- df_reg |> 
  mutate(
    est = paste0(
      round(estimate, 2), " (",
      round(conf.low, 2), ", ",
      round(conf.high, 2), ")"
    )
  ) |> 
  select(term, contrast, est, p.value, outcome, model, covariates) 

# Save data
# write.csv(df_reg, here("output", "table2b-regression_results.csv"), row.names = F)

# Effect modification by race and ethnicity ----
## Exclude multiracial Latino due to low reporting frequencies
df <- df |> filter(race_cat != "multiracial Latino")

## fit model for LAI over daily oral
### anticipated discrimination 
m1_a <- glm(prefer_inj_to_daily ~ indi_a_cat + race_cat + age + I(age^2) + educ_3cat + employed + income_3cat + sexor + unstable_house + prep_know_sum, data = df, family = binomial())
m2_a <- glm(prefer_inj_to_daily ~ indi_a_cat * race_cat + age + I(age^2) + educ_3cat + employed + income_3cat + sexor + unstable_house + prep_know_sum, data = df, family = binomial())

### daily discrimination 
m1_d <- glm(prefer_inj_to_daily ~ indi_d_cat + race_cat + age + I(age^2) + educ_3cat + employed + income_3cat + sexor + unstable_house + prep_know_sum, data = df, family = binomial())
m2_d <- glm(prefer_inj_to_daily ~ indi_d_cat * race_cat + age + I(age^2) + educ_3cat + employed + income_3cat + sexor + unstable_house + prep_know_sum, data = df, family = binomial())

### major discrimination
m1_m <- glm(prefer_inj_to_daily ~ indi_m_cat + race_cat + age + I(age^2) + educ_3cat + employed + income_3cat + sexor + unstable_house + prep_know_sum, data = df, family = binomial())
m2_m <- glm(prefer_inj_to_daily ~ indi_m_cat * race_cat + age + I(age^2) + educ_3cat + employed + income_3cat + sexor + unstable_house + prep_know_sum, data = df, family = binomial())

### likelihood ratio test
anova(m1_a, m2_a)
anova(m1_d, m2_d)
anova(m1_m, m2_m)

### combine results into a single table
### include robust standard errors
df_int1 <- bind_rows(
  avg_comparisons(m2_a, variables = "indi_a_cat", by = "race_cat", vcov = "HC3", comparison = "lnratioavg", transform = "exp"),
  avg_comparisons(m2_d, variables = "indi_d_cat", by = "race_cat", vcov = "HC3", comparison = "lnratioavg", transform = "exp"),
  avg_comparisons(m2_m, variables = "indi_m_cat", by = "race_cat", vcov = "HC3", comparison = "lnratioavg", transform = "exp")
) 

### format table
df_int1 <- df_int1 |> 
  mutate(
    est = paste0(
      round(estimate, 2), " (",
      round(conf.low, 2), ", ",
      round(conf.high, 2), ")"
    )
  ) |>
  select(term, contrast, race_cat, est) |>
  pivot_wider(names_from = race_cat, values_from = est)

## fit model for LAI over on-demand
### anticipated discrimination 
m1_a <- glm(prefer_inj_to_event ~ indi_a_cat + race_cat + age + I(age^2) + educ_3cat + employed + income_3cat + sexor + unstable_house + prep_know_sum, data = df, family = binomial())
m2_a <- glm(prefer_inj_to_event ~ indi_a_cat * race_cat + age + I(age^2) + educ_3cat + employed + income_3cat + sexor + unstable_house + prep_know_sum, data = df, family = binomial())

### daily discrimination 
m1_d <- glm(prefer_inj_to_event ~ indi_d_cat + race_cat + age + I(age^2) + educ_3cat + employed + income_3cat + sexor + unstable_house + prep_know_sum, data = df, family = binomial())
m2_d <- glm(prefer_inj_to_event ~ indi_d_cat * race_cat + age + I(age^2) + educ_3cat + employed + income_3cat + sexor + unstable_house + prep_know_sum, data = df, family = binomial())

### major discrimiantion 
m1_m <- glm(prefer_inj_to_event ~ indi_m_cat + race_cat + age + I(age^2) + educ_3cat + employed + income_3cat + sexor + unstable_house + prep_know_sum, data = df, family = binomial())
m2_m <- glm(prefer_inj_to_event ~ indi_m_cat * race_cat + age + I(age^2) + educ_3cat + employed + income_3cat + sexor + unstable_house + prep_know_sum, data = df, family = binomial())

### likelihood ratio test
anova(m1_a, m2_a)
anova(m1_d, m2_d)
anova(m1_m, m2_m)

### combine results into a single table
### include robust standard errors
df_int2 <- bind_rows(
  avg_comparisons(m2_a, variables = "indi_a_cat", by = "race_cat", vcov = "HC3", comparison = "lnratioavg", transform = "exp"),
  avg_comparisons(m2_d, variables = "indi_d_cat", by = "race_cat", vcov = "HC3", comparison = "lnratioavg", transform = "exp"),
  avg_comparisons(m2_m, variables = "indi_m_cat", by = "race_cat", vcov = "HC3", comparison = "lnratioavg", transform = "exp")
) 

# save results
# write.csv(df_int1, here("output", "table3-interaction1.csv"), row.names = F)
# write.csv(df_int2, here("output", "table3-interaction2.csv"), row.names = F)

### format table
df_int2 <- df_int2 |> 
  mutate(
    est = paste0(
      round(estimate, 2), " (",
      round(conf.low, 2), ", ",
      round(conf.high, 2), ")"
    )
  ) |>
  select(term, contrast, race_cat, est) |>
  pivot_wider(names_from = race_cat, values_from = est)

# Supplemental: prep knowledge items by race and ethnicity ----
v_names <- df |> 
  select(prep_know_1:prep_know_9) |> 
  names()

df_tab <- CreateTableOne(
  vars = v_names, 
  strata = "race_cat",
  data = df, 
  factorVars = v_names, 
  includeNA = T, 
  addOverall = T
)

tab <- print(
  df_tab, 
  quote = T,
  noSpaces = T, 
  printToggle = F
)

# Save table
# write.csv(tab, here("output", "sup_table1-prep_knowledge.csv"))
rm(v_names, df_tab, tab)