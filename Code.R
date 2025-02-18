library(devtools)  # Required for installing GitHub packages
# devtools::install_github("MathiasHarrer/dmetar") 

library(dmetar)    # Meta-analysis package
library(dplyr)     # Data manipulation
library(readr)     # Better CSV reading
library(meta)      # Meta-analysis computations

# Load dataset from GitHub
path <- "https://raw.githubusercontent.com/danialk20/Meta-analysis/refs/heads/main/Data.csv"
data <- read_csv(path)

# Convert categorical columns to factors
factor_cols <- c("Article_ID", "Instance_ID", "Objective", "Iterations")
data <- data %>% mutate(across(all_of(factor_cols), as.factor))

# Compute improvement percentage based on objective function direction
data <- data %>%
  mutate(Improvement = case_when(
    Objective == "MIN" ~ ((Without_SPSC - With_SPSC) / Without_SPSC) * 100,
    TRUE ~ ((With_SPSC - Without_SPSC) / Without_SPSC) * 100
  ))

# Compute study-level statistics
Studies <- data %>%
  group_by(Sk = Article_ID) %>%
  summarize(
    Ak = mean(Improvement, na.rm = TRUE),   # Mean improvement per study
    Sigmak = sd(Improvement, na.rm = TRUE), # Standard deviation
    Nk = n(),                               # Number of instances
    Vk = (Sigmak^2 / Nk),                   # Variance
    SESigmak = (Sigmak / sqrt(Nk))          # Standard error
  )

# Auxiliary calculations for Tau² estimation
Auxiliary <- Studies %>%
  summarize(
    sum_a = sum(1 / Vk, na.rm = TRUE),
    sum_b = sum(Ak / Vk, na.rm = TRUE),
    sum_c = sum((Ak^2) / Vk, na.rm = TRUE),
    sum_d = sum(1 / Vk^2, na.rm = TRUE)
  )

# Compute heterogeneity statistics
Df <- n_distinct(Studies$Sk) - 1  # Degrees of freedom
Q <- Auxiliary$sum_c - (Auxiliary$sum_b^2 / Auxiliary$sum_a)
C <- Auxiliary$sum_a - (Auxiliary$sum_d / Auxiliary$sum_a)
Tau2 <- max((Q - Df) / C, 0)  # Ensure Tau2 is non-negative

# Compute study weights
Studies <- Studies %>%
  mutate(
    Wk = 1 / (Vk + Tau2),        # Weight per study
    Wkn = Wk / sum(Wk),          # Normalized weight
    WAk = Wkn * Ak               # Weighted mean improvement
  )

# Perform meta-analysis using random effects model
m.gen_spsc <- metagen(
  TE = Ak, seTE = sqrt(Vk), studlab = Sk,
  data = Studies, sm = "SMD", fixed = FALSE,
  random = TRUE, method.tau = "DL", hakn = TRUE,
  title = "SP/SC"
)

# Display results
summary(m.gen_spsc)  # Meta-analysis summary
forest(m.gen_spsc)   # Forest plot
funnel(m.gen_spsc)   # Funnel plot
