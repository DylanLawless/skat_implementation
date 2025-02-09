# A genotype example ----
# install.packages('SKAT')
library(SKAT)
library(dplyr)
library(ggplot2)
# Set seed for reproducibility
set.seed(123)

# Generate random genetic data
n_subjects <- 200  # number of subjects
n_variants <- 10  # number of variants

# Simulate phenotypes (1: disease, 0: healthy)
phenotype <- rbinom(n_subjects, 1, 0.5)

# Simulate covariates: gender (binary), age (continuous), PC1 (continuous)
gender <- rbinom(n_subjects, 1, 0.5)
age <- rnorm(n_subjects, 50, 10)
PC1 <- rnorm(n_subjects, 0, 1)


# Simulate genotypes (binary: 0/1/2 copies of minor allele)
genotype <- matrix(rbinom(n_subjects * n_variants, 2, 0.5), nrow = n_subjects)

# Melt the genotype matrix for visualization
genotype_melted <- reshape2::melt(genotype)

# Change column names to "sample" and "variant"
colnames(genotype_melted)[1:2] <- c("sample", "variant")

# Plot the genotype matrix
library(scales)
p1 <- ggplot(genotype_melted, aes(y = sample, x = variant, fill = factor(value))) +
	geom_tile() +
	scale_fill_manual(values = c("grey90", "grey50", "grey0")) +
	theme_minimal() +
	labs(#x = "", y = "", 
		title = "Genotype matrix",
		fill = "Genotype\nvalue")

p1

# Simulate weights
weights <- c(rep(0.5, n_variants / 2), rep(2, n_variants / 2))

# Melt the genotype matrix for visualization
genotype_weights <- reshape2::melt(weights)
genotype_weights$variant <- rownames(genotype_weights) %>% as.numeric()

p2 <- ggplot(genotype_weights, aes(y = 0, x = variant, fill = factor(value))) +
	geom_tile() +
	scale_fill_manual(values = c("pink", "red")) +
	theme_minimal() +
	theme(axis.text.y = element_blank()) +
	labs(#x = "", y = "", 
		title = "Genotype weights",
		fill = "Genotype\nweights") +
	ylab("weight")

p2

# Weighted genotype matrix
GW <- genotype %*% diag(sqrt(weights))

# Melt the genotype matrix for visualization
GW_melted <- reshape2::melt(GW)

# Change column names to "sample" and "variant"
colnames(GW_melted)[1:2] <- c("sample", "variant")

# Plot the genotype matrix
library(scales)
p3 <- ggplot(GW_melted, aes(y = sample, x = variant, fill = value)) +
	geom_tile() +
	scale_fill_gradientn(colours = c("white", "blue", "red"), 
								values = rescale(c(0, 1, 2))) +
	theme_minimal() +
	
	labs(
		title = "Genotype weighted matrix",
		fill = "Genotype\nvalue\nweighted")
p3

library(gridExtra)
# Combine the plots using gridExtra
grid.arrange(p1, p2, p3, ncol = 1)

# Centering matrix
C <- diag(n_variants) - 1 / n_variants

# Melt the centering matrix for visualization
C_melted <- reshape2::melt(C)
colnames(C_melted)[1:2] <- c("variant_X", "variant_Y")

# Plot the centering matrix
p_c <- ggplot(C_melted, aes(x = variant_X, y = variant_Y, fill = value)) +
	geom_tile() +
	scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
								midpoint = 0, limit = c(-min(C),max(C)), space = "Lab", 
								name="Centering\nvalue") +
	theme_minimal() +
	labs(x = "", y = "", 
		  title = "Centering matrix")
p_c

# Compute kernel
K <- GW %*% C %*% t(GW)
# Compute kernel
K <- t(GW) %*% GW

# Melt the kernel matrix for visualization
K_melted <- reshape2::melt(K)
colnames(K_melted)[1:2] <- c("sample_X", "sample_Y")

# Plot the kernel matrix
p_k <- ggplot(K_melted, aes(x = sample_X, y = sample_Y, fill = value)) +
	geom_tile() +
	scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
								midpoint = 0, limit = c(-max(abs(K)),max(abs(K))), space = "Lab", 
								name="Kernel\nvalue") +
	theme_minimal() +
	labs(
		title = "Kernel matrix")

p_k

# Let's look at the first 5x5 block of the kernel matrix
print(K[1:5, 1:5])

# Melt the kernel matrix for visualization
K_melted <- reshape2::melt(K)

colnames(K_melted)[1:2] <- c("sample", "variant")

# Plot the kernel matrix
p_k <- ggplot(K_melted, aes(x = sample, y = variant, fill = value)) +
	geom_tile() +
	scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
								midpoint = 0, limit = c(-1,1), space = "Lab", 
								name="Kernel\nvalue") +
	theme_minimal() +
	theme(#axis.text.x = element_blank(),
		#axis.text.y = element_blank(),
		#axis.ticks = element_blank(),
		legend.position = "bottom") +
	labs(
		title = "Kernel matrix")

p_k

# Combine the plots
grid.arrange(p1, p2, p3, p_c, ncol = 1)

# Now show the example kernel
p_k

# Skat ----
# 1. Fit null model
null_model <- lm(phenotype ~ gender + age + PC1)

# 2. Create weighted genotype matrix
# compute weights based on minor allele frequency (maf)
weights <- ifelse(apply(genotype, 2, mean) > 0.5, 1-apply(genotype, 2, mean), apply(genotype, 2, mean))

# scale weights so that they sum up to the number of variants
weights <- weights / sum(weights) * n_variants

# create weighted genotype matrix
genotype_weighted <- sweep(genotype, 2, sqrt(weights), "*")

# 3. Compute kernel matrix
# Compute kernel
K <- GW %*% t(GW)

# 4. Perform variance component test
# extract residuals from null model
y <- null_model$residuals

# variance of y
Vy <- var(y)

# variance explained by genotypes
Vg <- crossprod(y, K %*% y) / (n_subjects - 1)


# under the null hypothesis Vg should be 0
# therefore, we compare Vg to 0 using a chi-squared test
p_value <- 1 - pchisq(Vg / Vy, df = 1)

print(paste("P-value:", p_value))

