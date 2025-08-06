library(dplyr)     # Data manipulation
library(readr)     # Better CSV reading
library(meta)      # Meta-analysis computations
library(stringr)   # Convert to list

# Load main dataset
path <- "https://raw.githubusercontent.com/danialk20/Meta-analysis/refs/heads/main/Data.csv"
data <- read_csv(path)

# Load Table 3 with additional information of the implementations
path_2 <- "https://raw.githubusercontent.com/danialk20/Meta-analysis/refs/heads/main/Table3.csv"
data_2 <- read_delim(path_2, delim = ",", show_col_types = FALSE)

# Merge characteristics to the Article_IDs
data <- data %>%
  inner_join(data_2, by = "Article_ID")

# Convert categorical columns to factors
factor_cols <- c("Article_ID", "Instance_ID", "Objective", "Iterations", "Iterative_Postopt", "Component", "Algorithm_Type")
data <- data %>% mutate(across(all_of(factor_cols), as.factor))

# Compute improvement percentage based on objective function direction
data <- data %>%
  mutate(Improvement = case_when(
    Objective == "MIN" ~ ((Without_SPSC - With_SPSC) / Without_SPSC) * 100,
    TRUE ~ ((With_SPSC - Without_SPSC) / Without_SPSC) * 100
  ))

# Compute study-level statistics
Studies <- data %>%
  group_by(Article_ID = Article_ID, Component = Component,
           Iterative_Postopt = Iterative_Postopt,
           Algorithm_Type = Algorithm_Type) %>%
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
Df <- n_distinct(Studies$Article_ID) - 1  # Degrees of freedom
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

# Random effects model
m.gen_spsc <- metagen(
  TE = Ak, seTE = sqrt(Vk), studlab = Article_ID,
  data = Studies, sm = "SMD", common = FALSE,
  random = TRUE, method.tau = "DL",
  title = "SP/SC"
)

# Display results
summary(m.gen_spsc)  # Meta-analysis summary
forest(m.gen_spsc)   # Forest plot
funnel(m.gen_spsc)   # Funnel plot


##################################################################
# Subgroup analysis - SP vs SC
##################################################################


# Random effects model
m.gen_sp_sc <- metagen(
  TE = Ak, seTE = sqrt(Vk), studlab = Article_ID,
  data = Studies, sm = "SMD", common = FALSE,
  random = TRUE, method.tau = "DL",
  title = "SP-SC", 
  subgroup = Component,
  tau.common = FALSE
)

# Display results
summary(m.gen_sp_sc)
forest(m.gen_sp_sc)


##################################################################
# Subgroup analysis - Iterative vs. Postoptimizer
##################################################################

# Random effects model
m.gen_it_po <- metagen(
  TE = Ak, seTE = sqrt(Vk), studlab = Article_ID,
  data = Studies, sm = "SMD", common = FALSE,
  random = TRUE, method.tau = "DL",
  title = "Iterative-Postoptimizer", 
  subgroup = Iterative_Postopt,
  tau.common = FALSE
)

# Display results
summary(m.gen_it_po)
forest(m.gen_it_po)


##################################################################
# Subgroup analysis - Algorithm Type
##################################################################

# Random effects model
m.gen_at <- metagen(
  TE = Ak, seTE = sqrt(Vk), studlab = Article_ID,
  data = Studies, sm = "SMD", common = FALSE,
  random = TRUE, method.tau = "DL",
  title = "Algorithm Type", 
  subgroup = Algorithm_Type,
  tau.common = FALSE
)

# Display results
summary(m.gen_at)
forest(m.gen_at)
