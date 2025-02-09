# Methods from other packages ----
# Requires functions:
# sqrt.inv, (not used here)
# saddle, (not used here)
# Liu.pval,  (not used here)
# Sadd.pval,  (not used here)
# aSKAT.b (not used here)
# KAT.pval
# aSKAT
# aSKAT.c
# davies

# Methods from SSKAT: Jun Chen (chen.jun2@mayo.edu) ----

# aSKAT
aSKAT <- function(formula.H0, data = NULL, K, type = c('binary', 'continuous'), acc = 0.00001, lim = 10000, tol = 1e-10) {
	type <- match.arg(type)
	if (type == 'continuous') {
		res <- aSKAT.c(formula.H0, data, K, acc = acc, lim = lim, tol = tol) 
	}
	if (type == 'binary') {
		res <- aSKAT.b(formula.H0, data, K, acc = acc, lim = lim, tol = tol) 
	}
	res
}

# sqrt.inv
sqrt.inv <- function (V2) {
	eig.obj <- eigen(V2, symmetric = TRUE)
	vectors <- eig.obj$vectors
	values <- eig.obj$values
	ind <- values >= 1e-10
	values <- values[ind]
	vectors <- vectors[, ind]
	temp <- t(vectors) / (values)
	Vi2 <- vectors  %*% temp
	temp <- t(vectors) / sqrt(values)
	Vi <- vectors  %*% temp
	return(list(Vi = Vi, Vi2 = Vi2, rank = length(values)))
}

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
}

# aSAKT.c
aSKAT.c <- function (formula.H0, data = NULL, K, acc = 0.00001, lim = 10000, tol = 1e-10) {
	res <- resid(lm(formula.H0, data))
	s2 <- sum(res^2) 
	X1 <- model.matrix(formula.H0, data)
	D0  <- diag(length(res))  
	P0  <- D0 - X1 %*% solve(t(X1) %*% X1) %*% t(X1)
	PKP <- P0 %*% K %*% P0
	q <- as.numeric(res %*% K %*% res / s2)
	ee <- eigen(PKP - q * P0, symmetric = T) 
	lambda <- ee$values[abs(ee$values) >= tol]
	p.value <- KAT.pval(0, lambda=sort(lambda, decreasing=T), acc = acc, lim = lim)
	return(list(p.value=p.value, Q.adj=q))
}

# Method: Davies from CompQuadForm----
library(CompQuadForm)
# https://github.com/cran/CompQuadForm
# However it required the qfc.cpp which is too long, and in C, to replicate.

# CUSTOM STANDALONE VERSION ----
# Custom method continuous ----

# Required libraries
library(CompQuadForm)

# Generate some data
set.seed(123)
Y <- rnorm(100) # phenotypic data
Z <- matrix(rnorm(200), 100, 2) # Covariates such as sex, age, PCs
G <- matrix(rbinom(1000, 1, 0.25), 100, 10) # Genotypic data
K <- G %*% t(G) # Kernel matrix quantifying similarities between samples
data <- data.frame(Y = Y, Z1 = Z[,1], Z2 = Z[,2]) # Dataframe containing the variables

# Linear regression formula to adjust variables
formula.H0 <- Y ~ Z1 + Z2

# Integration accuracy for Davies' method
acc = 0.00001

# Maximum number of integration terms for Davies' method
lim = 10000

# Eigenvalue cutoff, below which is considered to be 0, to reduce the computation burden
tol = 1e-10

# Residuals of the linear regression model
res <- resid(lm(formula.H0, data))

# Residual sum of squares
s2 <- sum(res^2) 

# Design matrix for the linear regression model
X1 <- model.matrix(formula.H0, data)

# Identity matrix
D0  <- diag(length(res))  

# Projection matrix for the null hypothesis
P0  <- D0 - X1 %*% solve(t(X1) %*% X1) %*% t(X1)

# Calculating the product of projection and kernel matrices
PKP <- P0 %*% K %*% P0

# Adjusted score statistic
q <- as.numeric(res %*% K %*% res / s2)

# Eigenvalues of PKP - q * P0
ee <- eigen(PKP - q * P0, symmetric = T) 

# Eigenvalues above the tolerance limit
lambda <- ee$values[abs(ee$values) >= tol]

# Calculating the p-value
p.value <- KAT.pval(0, lambda=sort(lambda, decreasing=T), acc = acc, lim = lim)

# Returning results as a list
res <- list(p.value = p.value, Q.adj = q)

# Print the result
print(res)



# The primary difference between SKAT and aSKAT lies in the sample size adjustment incorporated in aSKAT.
# 
# SKAT assumes that the sample size is sufficiently large such that the central limit theorem holds. Under the central limit theorem, the distribution of test statistics converges to a normal distribution as the sample size increases. Therefore, in large samples, SKAT can reliably use a chi-square distribution to approximate the null distribution of the test statistic.
# 
# However, when the sample size is small, the assumption of the central limit theorem may not hold, and the approximation of the null distribution using a chi-square distribution may not be accurate. This can result in inflated type I error rate (the probability of falsely rejecting the null hypothesis), leading to an excess of false-positive findings.
# 
# aSKAT addresses this issue by incorporating an adjustment for small sample sizes. Specifically, aSKAT computes the exact null distribution of the test statistic without relying on the central limit theorem, which allows it to accurately control the type I error rate regardless of the sample size. This makes aSKAT a more reliable method for studies with small sample sizes compared to SKAT.
# 
# Overall, the main conceptual difference between SKAT and aSKAT is how they handle small sample sizes. Both methods are designed to test for association between a set of genetic variants and a trait, but aSKAT provides a more robust approach when dealing with small sample sizes.

# The computation of the null distribution without relying on the central limit theorem occurs when the adjusted score statistic Q.adj is calculated and its associated p-value is derived using the Davies method.
# Here's a breakdown of the relevant steps in the script:
# Calculating the adjusted score statistic Q.adj:
# ``` r
# Q <- crossprod(resid.lm) # calculate the raw score statistic
# Q.adj <- sum(Q / eig.val) # adjust the score statistic based on the eigenvalues of the kernel matrix
# ```
# The adjusted score statistic Q.adj is calculated as the sum of the raw score statistic Q divided by the eigenvalues of the kernel matrix. This adjustment accounts for the correlation structure in the data, as captured by the kernel matrix.
# Calculating the p-value associated with Q.adj:
# ```r
# lambda <- (eig.val[eig.val > tol])^(-2) # calculate the eigenvalues for the Davies method
# p.value <- CompQuadForm::davies(Q.adj, lambda, acc = acc, lim = lim)$Qq # calculate the p-value
# ```
# The Davies method is used to calculate the p-value associated with Q.adj. This method computes the exact distribution of a quadratic form in normal variables, thereby not relying on the central limit theorem.
# Note that the Davies method is implemented in the CompQuadForm::davies function from the CompQuadForm package, and that package must be installed and loaded for this code to work.
# This exact computation of the null distribution allows aSKAT to accurately control the type I error rate regardless of the sample size, which is particularly important for studies with small sample sizes.
