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
	labs(
		title = "Genotype matrix",
		fill = "Genotype\nvalue"
	)

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
	labs(
		title = "Genotype weights",
		fill = "Genotype\nweights"
	) +
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
		fill = "Genotype\nvalue\nweighted"
	)

p3
# Centering matrix
C <- diag(n_variants) - 1 / n_variants

# Melt the centering matrix for visualization
C_melted <- reshape2::melt(C)
colnames(C_melted)[1:2] <- c("variant_X", "variant_Y")

# Plot the centering matrix
p_c <- ggplot(C_melted, aes(x = variant_X, y = variant_Y, fill = value)) +
	geom_tile() +
	scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
								midpoint = 0, limit = c(-min(C), max(C)), space = "Lab", 
								name="Centering\nvalue") +
	theme_minimal() +
	labs(
		title = "Centering matrix"
	)

p_c

# Compute kernel
K <- t(GW) %*% GW


# This computation creates a square symmetric matrix (the kernel), which represents pairwise similarities between samples (rows of the GW matrix). The multiplication by the transpose first (i.e., t(GW)) corresponds to calculating the dot product of each pair of samples, effectively comparing every individual to every other individual.

# Melt the kernel matrix for visualization
K_melted <- reshape2::melt(K)
colnames(K_melted)[1:2] <- c("sample_X", "sample_Y")

# Plot the kernel matrix
p_k <- ggplot(K_melted, aes(x = sample_X, y = sample_Y, fill = value)) +
	geom_tile() +
	scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
								midpoint = 0, limit = c(-max(abs(K)), max(abs(K))), space = "Lab", 
								name="Kernel\nvalue") +
	theme_minimal() +
	labs(
		title = "Kernel matrix"
	)

p_k
library(gridExtra)
# Combine the plots using gridExtra
grid.arrange(p1, p2, p3, p_c, p_k, ncol = 1)

# Skat ----
# 1. Fit null model
null_model <- glm(phenotype ~ gender + age + PC1, family=binomial)

# 2. Perform variance component test
# extract residuals from null model
y <- null_model$residuals

# variance of y
Vy <- var(y)

# For the computation of the genetic variance Vg, we should be using the kernel matrix K and the vector of residuals y from the null model. However, we need to be careful about the dimensions of these objects to ensure that the matrix multiplication is valid.
# Based on the typical SKAT method, Vg is calculated as the proportion of phenotypic variance that can be explained by the genetic variants (represented by the kernel K). A common method of estimating Vg in the context of a kernel-based association test is:

# Vg: variance explained by genotypes
Vg <- crossprod(y %*% GW) / (n_subjects - 1)

# This formula is computing the variance explained by the genotypes by calculating the sum of squared residuals when y is projected onto the space spanned by the genotype matrix GW. Then it is normalizing this by the degrees of freedom (n_subjects - 1) to get an estimate of the variance.





# Null model ----
my_SKAT_Null_Model = function(formula, data){
	
	# check missing 
	obj1 <- model.frame(formula, na.action = na.omit, data)
	obj2 <- model.frame(formula, na.action = na.pass, data)
	
	n <- dim(obj2)[1]
	n1 <- dim(obj1)[1]
	id_include <- complete.cases(obj1, obj2)
	
	if(n - n1 > 0){
		MSG <- sprintf("%d samples have missing phenotype or covariates. Excluded from the analysis!",n - n1)
		warning(MSG, call.=FALSE)
	}
	
	# Fit a linear model
	lm_fit <- lm(formula, data[id_include, ])
	
	re <- list(residuals = residuals(lm_fit), fitted.values = fitted.values(lm_fit), model = lm_fit)
	
	class(re) <- "my_SKAT_NULL_Model"
	re$n.all <- n
	return(re)
}





lm_fit <- glm(phenotype ~ gender + age + PC1, family=binomial)
re <- list(residuals = residuals(lm_fit), fitted.values = fitted.values(lm_fit), model = lm_fit)


The parameters in this function include:
	
res: The residuals from the null model.
Z: The matrix of genotypes, one row per individual and one column per variant.
X1: The matrix of covariates, one row per individual and one column per covariate.
kernel: The kernel to be used in the test.
weights: The weights for the genetic variants.
pi_1: The proportion of non-centrality parameters that are non-zero in the mixture model.
method: The method for combining the p-values.
res.out: The output from Get_ SKAT_Null_Model.
n.Resampling: The number of bootstrap resamplings to perform for the adjustment of small samples.
r.corr: The correlation parameter in the linear mixed model.
IsMeta: A logical value indicating whether to use Meta-SKAT or not.


SKAT.logistic.Linear = function(res,Z,X1, kernel, weights = NULL, pi_1, method,res.out,n.Resampling,r.corr, IsMeta=FALSE){
	
	
	if(length(r.corr) > 1 && dim(Z)[2] == 1){
		r.corr=0
	}
	
	if(IsMeta){
		
		# re = SKAT_RunFrom_MetaSKAT(res=res,Z=Z, X1=X1, kernel=kernel, weights=weights, pi_1=pi_1
		# 									, out_type="D", method=method, res.out=res.out, n.Resampling=n.Resampling, r.corr=r.corr)
		
	} else if(length(r.corr) == 1 ){
		
		# re = KMTest.logistic.Linear(res,Z,X1, kernel, weights, pi_1, method
		# 									 , res.out, n.Resampling, r.corr)
		
	} else {
		
		
		re =SKAT_Optimal_Logistic(res, Z, X1, kernel, weights, pi_1, method
										  , res.out, n.Resampling, r.corr)
		
	}
	
	return(re)
}





# Fit logistic regression model
lm_fit <- glm(phenotype ~ gender + age + PC1, family=binomial)

# Prepare the object for SKAT
re <- list(residuals = residuals(lm_fit), 
			  fitted.values = fitted.values(lm_fit), 
			  model = lm_fit)

# Rest of your variables
Z <- genotype
X1 <- cbind(gender, age, PC1)
kernel <- K
weights <- weights
pi_1 <- 0.5
method <- "optimal"
res.out <- lm_fit
n.Resampling <- 1000
r.corr <- 0
n <- 10
m <- 200
# Finally, call your function
result <- SKAT.logistic.Linear(res = re, Z = Z, X1 = X1, kernel = kernel, weights = weights, 
										 pi_1 = pi_1, method = method, res.out = res.out, 
										 n.Resampling = n.Resampling, r.corr = r.corr)

result


names(result)


# The SKAT.logistic.Linear function as you've provided it here does not actually perform the SKAT test, but rather it fits a logistic regression model to your data, depending on the method parameter it calls either SKAT_RunFrom_MetaSKAT, KMTest.logistic.Linear or SKAT_Optimal_Logistic.
# To perform the actual SKAT test and get a p-value, you will want to call the SKAT function from the SKAT package. Here is an example:

SKAT.linear.Linear = function(res,Z,X1, kernel, weights = NULL, s2, method,res.out,n.Resampling,r.corr, IsMeta=FALSE){
	
	if(length(r.corr) > 1 && dim(Z)[2] == 1){
		r.corr=0
	}
	
	if(IsMeta){
		re = SKAT_RunFrom_MetaSKAT(res=res,Z=Z, X1=X1, kernel=kernel, weights=weights, s2=s2
											, out_type="C", method=method, res.out=res.out, n.Resampling=n.Resampling, r.corr=r.corr)
		
		
	} else if(length(r.corr) == 1 ){
		
		re = KMTest.linear.Linear(res,Z,X1, kernel, weights, s2, method
										  , res.out, n.Resampling, r.corr)
		
	} else {
		
		
		
		re =SKAT_Optimal_Linear(res, Z, X1, kernel, weights, s2, method
										, res.out, n.Resampling, r.corr)
	}
	return(re)
}


KMTest.linear.Linear = function(res,Z,X1, kernel, weights, s2, method,res.out,n.Resampling,r.corr){
	
	
	n<-nrow(Z)
	# Weighted Linear Kernel 
	if (kernel == "linear.weighted") {
		Z = t(t(Z) * (weights))
	}
	
	if(r.corr == 1){
		Z<-cbind(rowSums(Z))
	} else if(r.corr > 0){
		
		p.m<-dim(Z)[2]	
		R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
		L<-chol(R.M,pivot=TRUE)
		Z<- Z %*% t(L) 
	}
	
	
	# get Q
	Q.Temp = t(res)%*%Z
	
	Q = Q.Temp %*% t(Q.Temp)/s2/2
	
	Q.res = NULL
	if(n.Resampling > 0){
		Q.Temp.res = t(res.out)%*%Z
		
		Q.res = rowSums(rbind(Q.Temp.res^2))/s2/2
	}
	
	W.1 = t(Z) %*% Z - (t(Z) %*%X1)%*%solve(t(X1)%*%X1)%*% (t(X1) %*% Z ) # t(Z) P0 Z
	
	
	if( method == "liu" ){
		
		out<-Get_Liu_PVal(Q, W.1, Q.res)    
		pval.zero.msg=NULL
		
	} else if( method == "liu.mod" ){
		
		out<-Get_Liu_PVal.MOD(Q, W.1, Q.res)    
		pval.zero.msg = NULL
		
	} else if( method == "davies"  ){
		
		out<-Get_Davies_PVal(Q, W.1, Q.res)    
		pval.zero.msg = out$pval.zero.msg
		
	} else {
		stop("Invalid Method!")
	}
	
	
	re<-list(p.value = out$p.value, p.value.resampling = out$p.value.resampling
				, Test.Type = method, Q = Q, param=out$param ,pval.zero.msg=pval.zero.msg )  
	
	return(re)
	
	
}





SKAT.linear = function(Z,y, X1, kernel = "linear", weights = NULL, method="liu", res.out=NULL, n.Resampling = 0, r.corr=r.corr){
	n = length(y) 
	m = ncol(Z) 
	mod = lm(y~X1 -1)
	s2 = summary(mod)$sigma**2
	res = mod$resid
	# If m >> p and ( linear or linear.weight) kernel than call 
	# Linear function
	if( (kernel =="linear" || kernel == "linear.weighted") && n > m){
		re = SKAT.linear.Linear(res,Z,X1, kernel, weights,s2,method,res.out,n.Resampling,r.corr)
	} else {  
		re = SKAT.linear.Other(res,Z,X1, kernel, weights,s2,method,res.out,n.Resampling)  
	}
	return(re)
}


# Your existing code...
lm_fit <- glm(phenotype ~ gender + age + PC1, family=binomial)
re <- list(residuals = residuals(lm_fit), fitted.values = fitted.values(lm_fit), model = lm_fit)

# Formulate the covariate matrix X1
X1 <- cbind(gender, age, PC1)

# Apply the SKAT.linear function
skat_linear_result <- SKAT.linear(Z = genotype, y = phenotype, X1 = X1, kernel = "linear", weights = weights, method = "liu", res.out = re, n.Resampling = 0, r.corr = NULL)

# skat_linear_result now contains the output from the SKAT.linear function
skat_linear_result























