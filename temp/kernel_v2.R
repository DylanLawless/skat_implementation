# # Load required packages
# library(kernlab)
# library(ggplot2)
# 
# # Sample dataset
# square_footage <- c(120, 250, 300, 150, 180, 240, 210, 200, 170, 280)  # Square footage of 10 houses
# price <- c(2000, 5200, 5900, 2500, 3100, 5100, 4500, 4200, 3500, 5800)  # Prices of 10 houses
# 
# # Data frame
# data <- data.frame(square_footage = square_footage, price = price)
# 
# # Order the dataset by square_footage for better plotting
# data <- data[order(data$square_footage),]
# 
# # Gaussian Kernel Regression
# gaussian_kernel_reg <- function(L) {
# 	# Perform the kernel regression
# 	model <- ksvm(price~., data = data, kernel = "rbfdot", kpar=list(sigma=L), type="eps-svr")
# 	
# 	# Predict on the training data
# 	predicted <- predict(model, data)
# 	
# 	# Plot
# 	ggplot() +
# 		geom_point(aes(x = square_footage, y = price), data = data, size = 3) +
# 		geom_line(aes(x = square_footage, y = predicted), data = data.frame(square_footage = data$square_footage, price = predicted), color = "red") +
# 		labs(title = paste("Gaussian Kernel, L =", L), x = "Square Footage (in hundreds)", y = "Price") +
# 		theme_minimal()
# }
# 
# # Laplace Kernel Regression
# laplace_kernel_reg <- function(L) {
# 	# Perform the kernel regression
# 	model <- ksvm(price~., data = data, kernel = "laplacedot", kpar=list(sigma=L), type="eps-svr")
# 	
# 	# Predict on the training data
# 	predicted <- predict(model, data)
# 	
# 	# Plot
# 	ggplot() +
# 		geom_point(aes(x = square_footage, y = price), data = data, size = 3) +
# 		geom_line(aes(x = square_footage, y = predicted), data = data.frame(square_footage = data$square_footage, price = predicted), color = "blue") +
# 		labs(title = paste("Laplace Kernel, L =", L), x = "Square Footage (in hundreds)", y = "Price") +
# 		theme_minimal()
# }
# 
# # Generate plots for different values of L
# gaussian_kernel_reg(0.05)
# gaussian_kernel_reg(1)
# gaussian_kernel_reg(20)
# 
# laplace_kernel_reg(0.05)
# laplace_kernel_reg(1)
# laplace_kernel_reg(20)



# Skat kernel ----

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


# The use of the kernel: In machine learning, kernels are used to measure the similarity or dissimilarity between each pair of data points. The kernel is a function that transforms the input data into a higher-dimensional space, which allows linear methods to be used for non-linear problems. In the case of SKAT, the kernel is used to capture the combined effects of multiple genetic variants.

# Computing the kernel: The kernel is computed by pre-multiplying the weighted genotype matrix GW with its transpose and the centering matrix C. The result is a 200x200 matrix because it reflects the pairwise similarity between all the subjects (200 in total) in the data. Each cell in the kernel matrix represents the "similarity" between two subjects' genotypes.

# Centering matrix
C <- diag(n_variants) - 1 / n_variants

# Melt the centering matrix for visualization
C_melted <- reshape2::melt(C)
colnames(C_melted)[1:2] <- c("variant_X", "variant_Y")

# Computing the Centering Matrix: Before computing the kernel, you create a centering matrix C. This is a special type of matrix used in multivariate statistics to mean-center an array of data, in this case, the genotype data. Mean-centering means subtracting the mean of the data from each data point, essentially centering the data around zero. This is important for many machine learning algorithms, including kernel methods.

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
p_k

# Simulating genetic data: You are generating binary genetic data for n_subjects and n_variants. You are also simulating the phenotype as well as a few covariates (gender, age, and PC1).
# Genotype matrix: You are visualizing the genotype matrix using ggplot2. Each cell in the heatmap represents a variant for a subject, colored based on the genotype value (0, 1, or 2).
# Simulating weights and Genotype weights visualization: You are defining weights and applying them to the variants. Half of the variants are given a weight of 0.5 and the other half a weight of 2. This is then visualized in a heatmap as well.
# Weighted Genotype matrix: You are creating a weighted genotype matrix by multiplying the genotype matrix with the square root of the weights. This weighted genotype matrix is then visualized in a heatmap as well.
# Kernel computation: You are computing the kernel matrix, which is essentially a measure of similarity between every pair of subjects based on their genotypes. The kernel is computed using the weighted genotype matrix and the centering matrix.
# Kernel matrix visualization: Finally, you are visualizing the kernel matrix. Each cell in this heatmap represents the similarity between a pair of subjects.


# Q. is this hypothetical kernel matrix then what is used in the SKAT method to separate genotypes that have different association with the outcome, such that kernel values ~ 1 would be one group, and kernel values ~ -1 would be another group, and this new information allows for the regression to be performed while correcting for the different direction of effects?



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
# K <- t(genotype_weighted) %*% genotype_weighted

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

