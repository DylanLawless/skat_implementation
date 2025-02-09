# A genotype example ----
# install.packages('SKAT')
library(SKAT)
library(dplyr)
library(ggplot2)

# Set up data ----
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

# Plot the genotype matrix ----
library(scales)
p_GM <- ggplot(genotype_melted, aes(y = sample, x = variant, fill = factor(value))) +
	geom_tile() +
	scale_fill_manual(values = c("grey90", "grey50", "grey0")) +
	theme_classic() +
	labs(
		title = "Genotype matrix",
		fill = "Genotype\nvalue"
	)

# p_GM
# Simulate weights
weights <- c(rep(0.5, n_variants / 2), rep(2, n_variants / 2))

# Melt the genotype matrix for visualization
genotype_weights <- reshape2::melt(weights)
genotype_weights$variant <- rownames(genotype_weights) %>% as.numeric()

# Plot the variant weights ----
p_W <- ggplot(genotype_weights, aes(y = 0, x = variant, fill = factor(value))) +
	geom_tile() +
	scale_fill_manual(values = c("pink", "red")) +
	theme_classic() +
	theme(axis.text.y = element_blank()) +
	labs(
		title = "Genotype weights",
		fill = "Genotype\nweights"
	) +
	ylab("weight")

# p_W
# Weighted genotype matrix
GW <- genotype %*% diag(sqrt(weights))

# Melt the genotype matrix for visualization
GW_melted <- reshape2::melt(GW)

# Change column names to "sample" and "variant"
colnames(GW_melted)[1:2] <- c("sample", "variant")

# Plot the weighted genotype matrix ----
library(scales)
p_GW <- ggplot(GW_melted, aes(y = sample, x = variant, fill = value)) +
	geom_tile() +
	scale_fill_gradientn(colours = c("white", "blue", "red"), 
								values = rescale(c(0, 1, 2))) +
	theme_classic() +
	
	labs(
		title = "Genotype weighted matrix",
		fill = "Genotype\nvalue\nweighted"
	)

# p_GW
# Centering matrix
C <- diag(n_variants) - 1 / n_variants

# Melt the centering matrix for visualization
C_melted <- reshape2::melt(C)
colnames(C_melted)[1:2] <- c("variant_X", "variant_Y")

# Plot the centering matrix ----
p_C <- ggplot(C_melted, aes(x = variant_X, y = variant_Y, fill = value)) +
	geom_tile() +
	scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
								midpoint = 0, limit = c(-min(C), max(C)), space = "Lab", 
								name="Centering\nvalue") +
	theme_classic() +
	labs(
		title = "Centering matrix"
	)

# p_c

library(gridExtra)
# Combine the plots using gridExtra
p_grid <- grid.arrange(p_GM, p_W, p_GW, p_C, ncol = 1)

ggsave(p_grid, filename = "./images/p_grid_GM_GW_GWM_CM.pdf", width = 5, height = 10)

# Compute kernel
# K <- t(GW) %*% GW

# Compute kernel
K <- GW %*% t(GW)

# Comment set : kernel

# Melt the kernel matrix for visualization
K_melted <- reshape2::melt(K)
colnames(K_melted)[1:2] <- c("sample_X", "sample_Y")

library(RColorBrewer)
# Define the YlOrRd color palette
ylorrd_colors <- brewer.pal(9, "YlOrRd")

# Plot the kernel matrix ----
p_k <- ggplot(K_melted, aes(x = sample_X, y = sample_Y, fill = value)) +
	geom_tile() +
	scale_fill_gradientn(colors = ylorrd_colors, name = "Kernel\nvalue") +
	theme_classic() +
	labs(
		title = "Kernel matrix"
	)

# p_k
ggsave(p_k, filename = "./images/p_kernel_1.pdf", width = 5, height = 4.5)

# Plot the kernel matrix clustered ----
# We could also cluster to see the features of the kernel matrix
# Modify the default color scheme
default_palette <- colorRampPalette(ylorrd_colors)
palette(default_palette)
heatmap(K, symm = TRUE)
heatmap_filepath <- "./images/p_kernel_2.pdf"
pdf(heatmap_filepath)
heatmap(K, symm = TRUE)
dev.off()

# Comment set: prompt 

# SKAT test ----
# Continuing from kernel_example.R with manual SKAT analysis
# Required library for p-value calculation
library(CompQuadForm)

# Phenotypic data
Y <- phenotype

# Covariates such as gender, age, PCs
Z <- data.frame(gender, age, PC1)
Z <- data.frame(gender, age)

# Create a data frame to contain the variables
data <- data.frame(Y = Y, Z)

# Linear regression model ----
# Create a linear regression formula to adjust variables
formula.H0 <- Y ~ gender + age + PC1
formula.H0 <- Y ~ gender + age
# Fit the linear regression model
linear_model <- lm(formula.H0, data)

# Calculate residuals
res <- resid(linear_model)

# Calculate residual sum of squares
s2 <- sum(res^2) 

# Calculate the design matrix for the linear regression model
X1 <- model.matrix(formula.H0, data)

# Identity matrix
D0 <- diag(length(res))  

# Create a dataframe from the matrix and melt it
D0_df <- reshape2::melt(D0)

#  Plot the identity matrix ----
p_D0 <- ggplot(D0_df, aes(x=Var1, y=Var2, fill=value)) + 
	geom_tile() +
	theme_classic() +
	labs(title = "Identity Matrix (D0)", x = "Row", y = "Column") +
	scale_fill_gradient(low = "white", high = "steelblue")

ggsave(p_D0, filename = "./images/p_D0.pdf", width = 5, height = 4.5)

# Calculate projection matrix for the null hypothesis
P0 <- D0 - X1 %*% solve(t(X1) %*% X1) %*% t(X1)

# Create a dataframe from the matrix and melt it
P0_df <- reshape2::melt(P0)

# Plot the projection matrix for the null hypothesis ----
p_P0 <- ggplot(P0_df, aes(x=Var1, y=Var2, fill=value)) + 
	geom_tile() +
	theme_classic() +
	labs(title = "Projection Matrix (P0)", x = "Row", y = "Column") +
	scale_fill_gradient(low = "white", high = "steelblue")

ggsave(p_P0, filename = "./images/p_P0.pdf", width = 5, height = 4.5)


# Plot the first set of values from the projection matrix for the null hypothesis to illustrate more clearly ----
p_P0_first <- P0_df %>% filter(Var2 == 1) %>%
	filter(Var1 > 1) %>%
	ggplot(aes(x=Var1, y=value, fill=value)) + 
	geom_bar(stat = "identity") +
	theme_classic() +
	labs(title = "Projection Matrix (P0)", x = "Row", y = "Column") +
	scale_fill_gradient(low = "black", high = "steelblue") 

ggsave(p_P0_first, filename = "./images/p_P0_first.pdf", width = 5, height = 4.5)

# Calculate the product of projection and kernel matrices
PKP <- P0 %*% K %*% P0

print(dim(P0))
print(dim(K))

# Create a dataframe from the matrix and melt it
PKP_df <- reshape2::melt(PKP)

# Plot the product of projection and kernel matrices ----
p_PKP <- ggplot(PKP_df, aes(x=Var1, y=Var2, fill=value)) + 
	geom_tile() +
	theme_classic() +
	labs(title = "Product of Projection and Kernel Matrices (PKP)", x = "Row", y = "Column") +
	scale_fill_gradient(low = "white", high = "steelblue")

ggsave(p_PKP, filename = "./images/p_PKP.pdf", width = 5, height = 4.5)

# Adjusted score statistic
q <- as.numeric(res %*% K %*% res / s2)

# Eigenvalues calculation ----
# Calculate eigenvalues of PKP - q * P0
ee <- eigen(PKP - q * P0, symmetric = TRUE) 

# Filter eigenvalues above the tolerance limit
tol = 1e-10
lambda <- ee$values[abs(ee$values) >= tol]

# Plot the lambda ----
p_lambda <- ggplot(data.frame(lambda = lambda, sample = 1:length(lambda)), aes(x = sample, y = lambda)) +
	geom_point(fill = "black", size = 1, alpha = 0.5) +
	theme_classic() +
	labs(title = "Contributions of Lambda to Test Statistic", x = "Sample", y = "Lambda")

ggsave(p_lambda, filename = "./images/p_lambda.pdf", width = 5, height = 4.5)


# P-value calculation ----
# Use Davies method for p-value calculation
acc = 0.00001
lim = 10000
# p_value <- CompQuadForm::davies(q, lambda, acc = acc, lim = lim)

# KAT.pval
#Compute the tail probability of 1-DF chi-square mixtures
KAT.pval <- function(Q.all, lambda, acc=1e-9,lim=1e6){
	pval = rep(0, length(Q.all))
	i1 = which(is.finite(Q.all))
	for(i in i1){
		tmp = davies(Q.all[i],lambda,acc=acc,lim=lim); pval[i] = tmp$Qq
		if((tmp$ifault>0)|(pval[i]<=0)|(pval[i]>=1)) pval[i] = Sadd.pval(Q.all[i],lambda)
	}
	return(pval)
} # KAT.pval used to compute the p-value, which internally uses davies() but also includes additional checks and adjustments.

# Calculating the p-value
p.value <- KAT.pval(0, lambda=sort(lambda, decreasing=T), acc = acc, lim = lim)

# Save and print the results
res <- list(p_value = p.value, Q_adj = q)
print(res)

# Comment set: about skat ----


# Comment set: lambda ----














# Comment set : kernel ----

# Note the following comment is not that important and the correct version is already in latex.

# The kernel matrix is used in the SKAT method to account for the fact that different variants may have effects in different directions. The purpose of the kernel matrix in this context is to create a similarity matrix that captures these relationships among individuals based on their genotypes.
# The SKAT method uses this kernel matrix to perform a variance-component test to examine the joint effect of multiple variants in a region (e.g., a gene) on a phenotype.
# Since the dimensions of K and P0 need to be the same for the multiplication to work, and both of these matrices are designed to capture relationships among individuals (not variants), it's important that K is calculated as a 200x200 matrix, as I've described in the previous response.
# Therefore, I suggest that we proceed with calculating K as GW %*% t(GW) and then perform P0 %*% K %*% P0.

# I am not 100% certain about K here. 

# Larson 2019: Instead of the 
# usual regression approach of explicitly defining basis
# functions to represent the effect of f (⋅), kernel methods
# approximate f (⋅) via an N × N kernel matrix K , such that
# Ki.j = κ (Gi, Gj) for some kernel function κ (⋅,⋅). A simple
# and popular kernel function is the linear‐weighted kernel,
# κ (Gi, Gj) = ∑M wk Gik Gjk, for some set of a priori k=1
# defined weights wk. 

# sample size N 
# continuous‐valued outcome y = (y ,..., y )′.
# variant set of interest of size M
# N × M matrix -> G 
# Simultaneously modeling the linear genetic effects for all M variants in G while adjusting for X may be conducted using ordinary multiple linear regression,

#  covariance between genetic markers k and l, for subjects i and j in the same

# The semiparametric regression approach, such that the covariates are modeled parametrically though the genetic markers are modeled in a nonparametric manner. The regression equation can be expressed as
# yi=α0 +Xi'α+f(Gi)+εi
# where f (⋅) is some unknown smooth function. 


# Comment set: prompt ----
# To continue with the SKAT example (part 2 of 2), we need to use the kernel (K) from the first script (1 of 2) kernel_example.R.
# We've generated data, phenotypes, genotypes, and covariates in the first script and constructed a kernel matrix. Now, we will run the SKAT using this kernel matrix and the function defined in skat_example.R.
# Let's continue from kernel_example.R with the SKAT analysis:

# Comment set: about skat ----
# In the context of the SKAT (Sequence Kernel Association Test), the Q statistic is the result of a quadratic form calculation involving the genotype data and the kernel matrix, K. This calculation generates a single summary statistic that represents the combined effect of all variants within a gene or region.

# In a nutshell, the Q statistic in this setting is a measure of the overall genetic association signal across a group of genetic variants in relation to the trait or disease under study.

# In the mathematical form, the Q statistic is derived from the formula Q = y'Ky, where y is a vector of individual genotype scores adjusted for covariates and K is the kernel matrix. The kernel matrix K represents the pairwise relationships between variants, which may account for their correlations or weights.

# When SKAT is performed, it essentially tests whether the Q statistic is larger than expected under the null hypothesis of no association, implying that the group of variants as a whole exhibits an association with the phenotype of interest. The p-value of the test is derived from the distribution of Q under the null hypothesis.

# Therefore, the Q statistic is a critical component of the SKAT methodology, as it captures the overall signal of association across multiple genetic variants. The p-value associated with this Q statistic is then used to determine statistical significance of the association.

# However, it's worth mentioning that interpretation of the Q statistic by itself can be complex, as its value depends on the scale of the input data, the choice of kernel, and the number and correlations of variants included in the test. Typically, the Q statistic is used in conjunction with the resulting p-value to make conclusions about genetic associations.


# Comment set: lambda ----
# It appears that the majority of the lambda values are positive, indicating a positive contribution to the test statistic. The larger positive values suggest stronger contributions from those specific eigenvalues.

# The negative value -4.20686 appears at indices 11 to 21. This indicates negative eigenvalues, which can arise due to noise or instability in the data or method used. 
