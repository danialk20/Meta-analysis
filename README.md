# Meta-analysis Repository

## Description

This repository contains data, R code, and results related to a meta-analysis study on optimization techniques for the Vehicle Routing Problem (VRP). The study focuses on the impact of using set partitioning/covering as a post-optimization step to enhance solutions obtained through meta- and matheuristics. The analysis was conducted using R and includes both forest and funnel plots to visualize the findings.

## Repository Structure

- **Data.csv**: The dataset in CSV format.
- **Data.xlsx**: The dataset in Excel format.
- **Code.R**: R script used for performing the meta-analysis.
- **Plots/**: Contains the generated forest and funnel plots.

## Meta-analysis Process

The meta-analysis was conducted using the `dmetar` and `metafor` packages in R. The steps involved:

1. **Loading Data**: The dataset includes multiple studies with performance metrics related to VRP optimization.
2. **Data Processing**: Converting necessary columns into factors and calculating improvement percentages.
3. **Statistical Computation**: Estimating effect sizes, standard errors, and between-study variance.
4. **Model Fitting**: Using a random-effects model to compute overall estimates.
5. **Visualization**: Generating forest and funnel plots for result interpretation.

### Background on the Meta-analysis

The Vehicle Routing Problem (VRP) is one of the most researched topics in Operations Research (OR). Over time, VRP research has shifted from solving the classical capacitated VRP to addressing Multi-Attribute VRPs (MAVRPs), which incorporate more realistic constraints. MAVRPs are NP-Hard problems, making exact optimization challenging. As a result, OR practitioners frequently employ meta- and matheuristic approaches to find high-quality solutions efficiently.

While many meta- and matheuristic approaches generate multiple feasible solutions, only a subset is typically considered for the final result. Some studies have explored leveraging all feasible solutions by decomposing them into individual routes and solving a set partitioning/covering problem to determine the best combination of routes. This approach can be applied as a post-optimization subroutine or iteratively during the execution of the metaheuristic.

This meta-analysis aims to quantify the effectiveness of incorporating set partitioning/covering techniques into VRP optimization. By analyzing existing studies, we assess whether these techniques lead to measurable improvements in solution quality.

## Results

- **Forest Plot**: Displays effect sizes and confidence intervals for each study.
- **Funnel Plot**: Assesses potential publication bias.

## Usage

To replicate the analysis:

1. Download and open `Code.R` in RStudio.
2. Ensure the required packages (`meta`, `dplyr`, `readr`) are installed.
3. Run the script to reproduce the results.

## Contact

For any questions or contributions, please reach out through GitHub.

