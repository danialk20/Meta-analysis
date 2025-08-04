library(dplyr)     # Data manipulation
library(readr)     # Better CSV reading
library(meta)      # Meta-analysis computations
library(stringr)   # Convert to list

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
  group_by(Article_ID = Article_ID) %>%
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

# Perform meta-analysis using random effects model
m.gen_spsc <- metagen(
  TE = Ak, seTE = sqrt(Vk), studlab = Article_ID,
  data = Studies, sm = "SMD", fixed = FALSE,
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

# Load Table 3 with methods and attributes
path_2 <- "https://raw.githubusercontent.com/danialk20/Meta-analysis/refs/heads/main/Table3.csv"
data_2 <- read_delim(path_2, delim = ",", show_col_types = FALSE)

# Subset original data into SP and SC components
data_SP <- data %>%
  inner_join(data_2, by = "Article_ID") %>%
  filter(Component == "SP")

data_SC <- data %>%
  inner_join(data_2, by = "Article_ID") %>%
  filter(Component == "SC")

###############################################
# Meta-analysis

# Compute study-level statistics for SP
Studies_SP <- data_SP %>%
  group_by(Article_ID = Article_ID) %>%
  summarize(
    Ak = mean(Improvement, na.rm = TRUE),
    Sigmak = sd(Improvement, na.rm = TRUE),
    Nk = n(),
    Vk = (Sigmak^2 / Nk),
    SESigmak = (Sigmak / sqrt(Nk))
  )

# Meta-analysis for SP
m.gen_sp <- metagen(
  TE = Ak, seTE = sqrt(Vk), studlab = Article_ID,
  data = Studies_SP, sm = "SMD", fixed = FALSE,
  random = TRUE, method.tau = "DL",
  title = "SP/SC"
)

# Output results for SP
summary(m.gen_sp)
forest(m.gen_sp)

# Compute study-level statistics for SC
Studies_SC <- data_SC %>%
  group_by(Article_ID = Article_ID) %>%
  summarize(
    Ak = mean(Improvement, na.rm = TRUE),
    Sigmak = sd(Improvement, na.rm = TRUE),
    Nk = n(),
    Vk = (Sigmak^2 / Nk),
    SESigmak = (Sigmak / sqrt(Nk))
  )

# Meta-analysis for SC
m.gen_sc <- metagen(
  TE = Ak, seTE = sqrt(Vk), studlab = Article_ID,
  data = Studies_SC, sm = "SMD", fixed = FALSE,
  random = TRUE, method.tau = "DL",
  title = "SP/SC"
)

# Output results for SC
summary(m.gen_sc)
forest(m.gen_sc)

###############################################
# Difference in means

# Welch Two Sample t-test
t_test_component <- t.test(
  data_SP$Improvement,
  data_SC$Improvement,
  var.equal = FALSE,
  alternative = "two.sided"
)

print(t_test_component)


##################################################################
# Subgroup analysis - Attributes
##################################################################

# Split string columns into lists
data_2 <- data_2 %>%
  mutate(
    Method = strsplit(Method, "~"),
    Attributes = strsplit(Attributes, "~")
  )

# Define helper function to generate subgroup summaries
generate_subgroup_summary <- function(column_name, tag) {
  unique_values <- unique(unlist(data_2[[column_name]]))
  results <- data.frame()
  
  for (val in unique_values) {
    df_tmp <- data_2[sapply(data_2[[column_name]], function(x) val %in% x), ] %>%
      inner_join(data, by = "Article_ID") %>%
      select(Article_ID, Improvement) %>%
      mutate(Group = val)
    
    results <- bind_rows(results, df_tmp)
  }
  
  results %>%
    group_by(Group) %>%
    summarize(
      M = mean(Improvement, na.rm = TRUE),
      SD= sd(Improvement, na.rm = TRUE),
      N = n(),
      .groups = "drop"
    ) %>%
    rename(!!tag := Group)
}

# Compute results by Attribute
results_attributes <- generate_subgroup_summary("Attributes", "Attribute")

print(results_attributes)


##################################################################
# Subgroup analysis - Methods
##################################################################

# Compute results by Method
results_methods <- generate_subgroup_summary("Method", "Method")

print(results_methods)


##################################################################
# Subgroup analysis - Algorithm Type
##################################################################

# Subset original data into Constructive and Local search
data_C <- data %>%
  inner_join(data_2, by = "Article_ID") %>%
  filter(Algorithm_Type == "Constructive")

data_LS <- data %>%
  inner_join(data_2, by = "Article_ID") %>%
  filter(Algorithm_Type == "Local search")

# Compute study-level statistics for Constructive
Studies_C <- data_C %>%
  group_by(Article_ID = Article_ID) %>%
  summarize(
    Ak = mean(Improvement, na.rm = TRUE),
    Sigmak = sd(Improvement, na.rm = TRUE),
    Nk = n(),
    Vk = (Sigmak^2 / Nk),
    SESigmak = (Sigmak / sqrt(Nk))
  )

# Meta-analysis for Constructive
m.gen_c <- metagen(
  TE = Ak, seTE = sqrt(Vk), studlab = Article_ID,
  data = Studies_C, sm = "SMD", fixed = FALSE,
  random = TRUE, method.tau = "DL",
  title = "SP/SC"
)

# Output results for Constructive
summary(m.gen_c)
forest(m.gen_c)

# Compute study-level statistics for Local Search
Studies_LS <- data_LS %>%
  group_by(Article_ID = Article_ID) %>%
  summarize(
    Ak = mean(Improvement, na.rm = TRUE),
    Sigmak = sd(Improvement, na.rm = TRUE),
    Nk = n(),
    Vk = (Sigmak^2 / Nk),
    SESigmak = (Sigmak / sqrt(Nk))
  )

# Meta-analysis for Local Search
m.gen_ls <- metagen(
  TE = Ak, seTE = sqrt(Vk), studlab = Article_ID,
  data = Studies_LS, sm = "SMD", fixed = FALSE,
  random = TRUE, method.tau = "DL",
  title = "SP/SC"
)

# Output results for LS
summary(m.gen_ls)
forest(m.gen_ls)

###############################################
# Difference in means

# Welch Two Sample t-test
t_test_algorithm <- t.test( 
  data_C$Improvement,
  data_LS$Improvement,
  var.equal = FALSE,
  alternative = "two.sided"
)

print(t_test_algorithm)