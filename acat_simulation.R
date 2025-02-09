# Loading necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(scales)

# Set seed for reproducibility
set.seed(123)

# Number of variants
n_variants <- 50

# Simulate Associated Variables
CADD_score <- runif(n_variants, 1, 20)  # Score can be any value
allele_frequency <- runif(n_variants, 0, 1)  # Score can be any value
pathogenicity <- rbinom(n_variants, 1, 0.5)  # Score can be any value

# Simulate P-values with a few outliers
generate_pvalues_with_outliers <- function(n, outlier_indices) {
  p_values <- runif(n, 0, 1)
  # Assigning very small values to selected outliers
  p_values[outlier_indices] <- runif(length(outlier_indices), 0, 0.001)
  return(p_values)
}

# Indices of potential outliers
outlier_indices <- sample(1:n_variants, 2)  # Randomly select 2 variants to be outliers

# Generate P-values
CADD_score_pval <- generate_pvalues_with_outliers(n_variants, outlier_indices)
allele_frequency_pval <- generate_pvalues_with_outliers(n_variants, outlier_indices)
pathogenicity_pval <- generate_pvalues_with_outliers(n_variants, outlier_indices)



# Combine variables into a single dataframe for visualization
variant_data <- data.frame(
  variant = rep(1:n_variants, times = 6),
  measure = rep(c("CADD_score", "allele_frequency", "pathogenicity"), each = n_variants * 2),
  type = rep(c(rep("Original Score", n_variants), rep("P-Value", n_variants)), times = 3),
  value = c(CADD_score, CADD_score_pval, allele_frequency, allele_frequency_pval, pathogenicity, pathogenicity_pval)  # Alternating Original Scores and P-values for each measure
)


# Melt the data for ggplot
variant_data_melted <- melt(variant_data, id.vars = c('variant', 'measure', 'type'))

# Calculate the significance threshold for P-values (not log-transformed)
significance_threshold <- 0.05 / n_variants

library(gridExtra)

# Separate the data into two parts: one for original scores and one for P-values
original_scores_data <- variant_data_melted %>% filter(type == "Original Score")
p_values_data <- variant_data_melted %>% filter(type == "P-Value")

# Plot for original scores
p_original_scores <- ggplot(original_scores_data, aes(x = variant, y = value, color = measure)) +
  geom_point() +
  theme_classic() +
  labs(
    title = "Original Scores",
    y = "Values",
    color = "Variables"
  )

# Plot for P-values with log scale
p_p_values <- ggplot(p_values_data, aes(x = variant, y = -log10(value), color = measure)) +
  geom_point() +
  theme_classic() +
  labs(
    title = "P-values",
    y = "-log10(P-values)"
  ) +
  geom_hline(aes(yintercept = -log10(significance_threshold)), 
             linetype = "dashed", color = "red")

# Combine the two plots
p_combined <- grid.arrange(p_original_scores, p_p_values, ncol = 1)
print(p_combined)

# Save the combined plot (optional)
ggsave("./images/p_combined.pdf", plot = p_combined) #, width = 8, height = 12)

# ACAT library ----
# The Aggregated Cauchy Association Test (ACAT) described in the paper is a method for combining individual P-values to assess the overall significance of a set of statistical tests. This approach is particularly relevant in large-scale data analyses, such as genome-wide association studies (GWAS), where you're testing many hypotheses simultaneously. Let's break down the method and the paper's main points:
# 
# ### Simplified Explanation
# 
# 1. **Combining P-values**: When you perform many statistical tests, you get a P-value for each test. ACAT combines these P-values into one overall P-value. This helps to understand if there's a significant effect when looking at all tests together.
# 
# 2. **Use of the Cauchy Distribution**: The method uses the Cauchy distribution to transform individual P-values. This transformation is chosen because it's mathematically convenient and works well under various conditions, including when tests are correlated or when only a few tests are significant (sparsity).
# 
# 3. **Calculation Simplicity**: The ACAT method results in a simple test statistic that's easy to calculate and interpret. This is crucial for handling massive datasets, as complex calculations can be computationally expensive.
# 
# 4. **Non-Asymptotic Theory**: The paper provides a theoretical foundation showing that the distribution of the ACAT test statistic under the null hypothesis (i.e., when there’s no real effect) can be approximated by a Cauchy distribution. This approximation holds true even under various dependency structures among tests.
# 
# 5. **Power and Accuracy**: The method is shown to be powerful in detecting true effects (especially under sparsity) and accurate in its P-value calculations. This means it's effective in identifying significant results and controlling false positives.
# 
# 6. **Application to GWAS**: The method is particularly suited for GWAS, where many genetic variants are tested for association with diseases. The paper demonstrates its use in a study of Crohn's disease.
# 
# ### Technical Understanding
# 
# 1. **ACAT Function**: The `ACAT` function you've provided in R takes a vector or matrix of P-values and optionally weights for each P-value. If no weights are given, it assumes equal weights.
# 
# 2. **Transforming P-values**: Each P-value is transformed using the Cauchy distribution. The transformation is `tan((0.5 - Pvals) * pi)`. This transformation is key to the method's effectiveness.
# 
# 3. **Weighted Sum**: The transformed P-values are then combined in a weighted sum. The weights can be customized based on the importance or reliability of each test.
# 
# 4. **P-value Calculation**: The resulting sum is then used to calculate the overall P-value of the test using the Cauchy distribution. This is a straightforward process compared to many other combination methods.
# 
# 5. **ACAT-V for Genetic Data**: `ACAT_V` is a specialized version of ACAT for genetic data, where you're combining P-values from genetic variant tests. It accounts for factors like minor allele count (MAC) and uses either the Burden test or the Cauchy method based on MAC thresholds.
# 
# In summary, ACAT is a powerful and computationally efficient method for combining multiple P-values, well-suited for large-scale datasets where traditional methods might struggle. Its applicability to GWAS and its theoretical robustness make it a valuable tool in statistical genetics and other fields involving massive data analysis.

# library(devtools)
# devtools::install_github("yaowuliu/ACAT")
# library(ACAT)
# 
# # ACAT Aggregated Cauchy Assocaition Test ----
# # Description
# A p-value combination method using the Cauchy distribution.
# 
# # Usage
# ACAT(Pvals, weights = NULL, is.check = TRUE)
# 
# # Arguments
# Pvals: a numeric vector/matrix of p-values. When it is a matrix, each column of p- values is combined by ACAT.
# weights: a numeric vector/matrix of non-negative weights for the combined p-values. When it is NULL, the equal weights are used.
# is.check: logical. Should the validity of Pvals (and weights) be checked? When the size of Pvals is large and one knows Pvals is valid, then the checking part can be skipped to save memory.
# 
# # Value: The p-value(s) of ACAT.
# 
# # Exmples
# p.values<-c(2e-02,4e-04,0.2,0.1,0.8)
# ACAT(Pvals=p.values)
# ACAT(matrix(runif(1000),ncol=10))
# 
# 
# # ACAT_V A set-based test that uses ACAT to combine the variant-level p-values. ----
# # Description
# A set-based test that uses ACAT to combine the variant-level p-values.
# 
# # Usage
# ACAT_V(G, obj, weights.beta = c(1, 25), weights = NULL, mac.thresh = 10)
# 
# # Arguments
# G: a numeric matrix or dgCMatrix with each row as a different individual and each column as a separate gene/snp. Each genotype should be coded as 0, 1, 2.
# obj: an output object of the NULL_Model function.
# weights.beta: a numeric vector of parameters for the beta weights for the weighted kernels. If you want to use your own weights, please use the “weights” parameter. It will be ignored if “weights” parameter is not null.
# weights: a numeric vector of weights for the SNP p-values. When it is NULL, the beta weight with the “weights.beta” parameter is used.
# mac.thresh: a threshold of the minor allele count (MAC). The Burden test will be used to aggregate the SNPs with MAC less than this thrshold.
# 
# # Details
# The Burden test is first used to aggregate very rare variants with Minor Allele Count (MAC) < mac.thresh (e.g., 10), and a Burden p-value is obtained. For each of the variants with MAC >= mac.thresh, a variant-level p-value is calculated. Then, ACAT is used to combine the variant-level p-values and the Burden test p-value of very rare variants.
# If weights.beta is used, then the weight for the Burden test p-value is demetermined by the average Minor Allele Frequency (MAF) of the variants with MAC < mac.thresh; if the user-specified weights is used, then the weight for the Burden test p-value is the average of weights of the variants with MAC < mac.thresh.
# Note that the weights here are for the SNP p-vlaues. In SKAT, the weights are for the SNP score test statistics. To transfrom the SKAT weights to the weights here, one can use the formula that weights = (skat_weights)^2*MAF*(1-MAF).
# 
# # Value
# The p-value of ACAT-V.
# 
# # Examples
# library(Matrix)
# data(Geno)
# G<-Geno[,1:100] # Geno is a dgCMatrix of genotypes
# Y<-rnorm(nrow(G))
# Z<-matrix(rnorm(nrow(G)*4),ncol=4)
# obj<-NULL_Model(Y,Z)
# ACAT_V(G,obj)
# 
# # NULL_Model Get parameters and residuals from the NULL model ----
# # Description
# Compute model parameters and residuals for ACAT-V
# 
# # Usage
# NULL_Model(Y, Z = NULL)
# Arguments
# Y: a numeric vector of outcome phenotypes.
# Z: a numeric matrix of covariates. Z must be full-rank. Do not include intercept in Z. The intercept will be added automatically.
# 
# # Details
# Y could only be continuous or binary. If Y is continuous, a linear regression model is fitted. If Y is binary, it must be coded as 0,1 and a logistic model is fitted.
# 
# # Value
# This function returns an object that has model parameters and residuals of the NULL model of no association between genetic variables and outcome phenotypes. After obtaining it, please use ACAT_V function to conduct the association test.
# 
# # Examples
# Y<-rnorm(10000)
# Z<-matrix(rnorm(10000*4),ncol=4)
# obj<-NULL_Model(Y,Z)

# work ----
# 
# Assuming variant_data_melted is the dataframe created in your previous code
# Filter the data to include only P-values
p_values_data <- variant_data_melted %>% 
  filter(type == "P-Value") %>% 
  mutate(log_p_value = -log10(value))

# Plot the histogram of P-values
p_histogram <- ggplot(p_values_data, aes(x = value)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Histogram of P-values for Variants",
    x = "P-values",
    y = "Frequency"
  )

# Display the histogram
print(p_histogram)

# QQ-Plot of P-values
p_qqplot <- ggplot(p_values_data, aes(sample = value)) +
  geom_qq() +
  geom_qq_line() +
  theme_minimal() +
  labs(
    title = "QQ-Plot of P-values for Variants",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )

# Display the QQ-plot
print(p_qqplot)

# Combine the two plots
p_hist_qq <- grid.arrange(p_histogram, p_qqplot, ncol = 1)
print(p_hist_qq)
ggsave("./images/p_hist_qq.pdf", plot = p_hist_qq) 

# install.packages("ACAT")  # Uncomment if ACAT is available on CRAN or another repository

library(ACAT)

# Assuming p_values_data is already filtered to include only P-values
p_values_per_variant <- p_values_data %>%
  group_by(variant) %>%
  summarize(aggregated_pval = mean(value))  # Aggregating P-values for each variant

# Apply ACAT to aggregated P-values for each variant
acat_results_per_variant <- apply(p_values_per_variant[, -1], 1, function(pvals) ACAT(pvals))

# Prepare data for visualization
acat_results_df <- data.frame(
  Variant = p_values_per_variant$variant,
  ACAT_P_Value = acat_results_per_variant
)

# Plot ACAT results for each variant
p_acat <- ggplot(acat_results_df, aes(x = Variant, y = -log10(ACAT_P_Value))) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "ACAT P-values for Each Variant",
    y = "-log10(ACAT P-value)",
    x = "Variant"
  )  +
  geom_hline(aes(yintercept = -log10(significance_threshold)), 
             linetype = "dashed", color = "red")

p_acat

ggsave("./images/p_acat.pdf", plot = p_acat) 

