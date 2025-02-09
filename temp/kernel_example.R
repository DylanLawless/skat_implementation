library(ggplot2)
library(plotly)
library(caret)

feature_map_1 <- function(X) {
	return(data.frame(X[,1], X[,2], X[,1]^2 + X[,2]^2))
}

library(MASS)

set.seed(123)
data <- data.frame(mvrnorm(n = 100, mu = c(0,0), Sigma = matrix(c(1,0,0,1), nrow = 2)))
data$y <- ifelse((data$X1^2 + data$X2^2) < 0.5, "Class 1", "Class 2")

Z <- feature_map_1(data[,1:2])

# Add Z to data
data$Z1 <- Z[,1]
data$Z2 <- Z[,2]
data$Z3 <- Z[,3]

# 2D scatter plot
ggplot(data, aes(x = X1, y = X2, color = y)) +
	geom_point() +
	labs(x = expression(x[1]), y = expression(x[2]), title = "Original dataset") +
	theme_minimal()

# 3D scatter plot

plot_ly((data |> filter(Z3 < 8)), x = ~Z1, y = ~Z2, z = ~Z3, color = ~y, colors = "viridis") %>%
	add_markers() %>%
	layout(scene = list(xaxis = list(title = expression(z[1])),
							  yaxis = list(title = expression(z[2])),
							  zaxis = list(title = expression(z[3]))),
			 title = "Transformed dataset")

# Note that we define a new feature map feature_map_2 that maps each 2-dimensional point to a 3-dimensional point, where the third dimension is computed as exp(-(X[,1]^2 + X[,2]^2)). The rest of the code is similar to the previous example.




# version 2 ----
# install.packages("e1071")

library(ggplot2)
library(plotly)
library(caret)
library(MASS)
library(e1071)

feature_map_1 <- function(X) {
	return(data.frame(X[,1], X[,2], X[,1]^2 + X[,2]^2))
}

set.seed(123)
data <- data.frame(mvrnorm(n = 100, mu = c(0,0), Sigma = matrix(c(1,0,0,1), nrow = 2)))
data$y <- ifelse((data$X1^2 + data$X2^2) < 0.5, "Class 1", "Class 2")

Z <- feature_map_1(data[,1:2])

# Add Z to data
data$Z1 <- Z[,1]
data$Z2 <- Z[,2]
data$Z3 <- Z[,3]

# 2D scatter plot
ggplot(data, aes(x = X1, y = X2, color = y)) +
	geom_point() +
	labs(x = expression(x[1]), y = expression(x[2]), title = "Original dataset") +
	theme_minimal()

# 3D scatter plot
plot_ly((data |> filter(Z3 < 8)), x = ~Z1, y = ~Z2, z = ~Z3, color = ~y, colors = "viridis") %>%
	add_markers() %>%
	layout(scene = list(xaxis = list(title = expression(z[1])),
							  yaxis = list(title = expression(z[2])),
							  zaxis = list(title = expression(z[3]))),
			 title = "Transformed dataset")

# Convert y to a factor with numeric levels
data$y <- as.factor(ifelse(data$y == "Class 1", 1, 2))

# SVM using linear kernel
model <- svm(y ~ Z1 + Z2 + Z3, data = data, kernel = "linear", type = "C-classification")

# Print model coefficients
print(model$coefs)





library(kernlab)
library(ggplot2)

# Define feature map
feature_map_3 <- function(X) {
	data.frame(sqrt(2) * X[,1] * X[,2], X[,1]^2, X[,2]^2)
}

# Generate dataset
set.seed(123)
X <- matrix(rnorm(200), ncol = 2)
y <- ifelse(X[,1]^2 + X[,2]^2 < 0.5, 1, -1)

# Apply feature map
Z <- feature_map_3(X)

# Define custom kernel function
my_kernel <- function(a, b) {
	Z_a <- feature_map_3(a)
	Z_b <- feature_map_3(b)
	tcrossprod(Z_a, Z_b)
}


# Create SVM with polynomial kernel of degree 2
model <- ksvm(y ~ ., data = data.frame(X), type = "C-svc", kernel = "polydot", kpar = list(degree = 2, scale = 1, offset = 0))

# Print accuracy score
print(mean(predict(model, X) == y))



# Plot decision boundary
ggplot(data.frame(X), aes(X1, X2)) +
	geom_point(aes(color = factor(y))) +
	 stat_contour(aes(z = predict(model, X)), geom = "polygon", bins = 1, alpha = 0.3) +
	labs(x = expression(x[1]), y = expression(x[2]), title = "Support Vector Machine with custom kernel") +
	theme_minimal()

x













# Figure 2 ----
library(ggplot2)
library(plotly)

feature_map_2 <- function(X) {
	return(data.frame(X[,1], X[,2], exp(-(X[,1]^2 + X[,2]^2))))
}

library(MASS)

set.seed(123)
data <- data.frame(mvrnorm(n = 100, mu = c(0,0), Sigma = matrix(c(1,0,0,1), nrow = 2)))
data$y <- ifelse((data$X1^2 + data$X2^2) < 0.5, "Class 1", "Class 2")

Z <- feature_map_2(data[,1:2])

# Add Z to data
data$Z1 <- Z[,1]
data$Z2 <- Z[,2]
data$Z3 <- Z[,3]

# 2D scatter plot
ggplot(data, aes(x = X1, y = X2, color = y)) +
	geom_point() +
	labs(x = expression(x[1]), y = expression(x[2]), title = "Original dataset") +
	theme_minimal()

# 3D scatter plot
plot_ly(data, x = ~Z1, y = ~Z2, z = ~Z3, color = ~y, colors = "viridis") %>%
	add_markers() %>%
	layout(scene = list(xaxis = list(title = expression(z[1])),
							  yaxis = list(title = expression(z[2])),
							  zaxis = list(title = expression(z[3]))),
			 title = "Transformed dataset")



# Figure 3 ----

library(plotly)
library(reshape2)

#load data

my_df <- iris
petal_lm <- lm(Petal.Length ~ 0 + Sepal.Length + Sepal.Width,data = my_df)

#Graph Resolution (more important for more complex shapes)
graph_reso <- 0.05

#Setup Axis
axis_x <- seq(min(my_df$Sepal.Length), max(my_df$Sepal.Length), by = graph_reso)
axis_y <- seq(min(my_df$Sepal.Width), max(my_df$Sepal.Width), by = graph_reso)

#Sample points
petal_lm_surface <- expand.grid(Sepal.Length = axis_x,Sepal.Width = axis_y,KEEP.OUT.ATTRS = F)
petal_lm_surface$Petal.Length <- predict.lm(petal_lm, newdata = petal_lm_surface)
petal_lm_surface <- acast(petal_lm_surface, Sepal.Width ~ Sepal.Length, value.var = "Petal.Length") #y ~ x

hcolors=c("red","blue","green")[my_df$Species]
iris_plot <- plot_ly(my_df, 
							x = ~Sepal.Length, 
							y = ~Sepal.Width, 
							z = ~Petal.Length,
							text = ~Species, # EDIT: ~ added
							type = "scatter3d", 
							mode = "markers",
							marker = list(color = hcolors))

	
	iris_plot <- add_trace(p = iris_plot,
								  z = petal_lm_surface,
								  x = axis_x,
								  y = axis_y,
								  type = "surface")

iris_plot



# SVD ----
library(e1071)
library(ggplot2)
library(plotly)

# Function to map the features to a higher-dimensional space
feature_map <- function(X) {
	return(data.frame(X[,1], X[,2], X[,1]^2 + X[,2]^2))
}

# Set seed for reproducibility
set.seed(123)

# Generate the data
data <- data.frame(mvrnorm(n = 100, mu = c(0,0), Sigma = matrix(c(1,0,0,1), nrow = 2)))
data$y <- ifelse((data$X1^2 + data$X2^2) < 0.5, 1, -1)

# Map the features
Z <- feature_map(data[,1:2])

# Add Z to data
data$Z1 <- Z[,1]
data$Z2 <- Z[,2]
data$Z3 <- Z[,3]

data$y <- as.factor(ifelse(data$y == "1", 1, 2))

# Perform the SVM classification
svm_model <- svm(y ~ ., data = data, kernel = "linear", cost = 10, scale = FALSE, class.weights = c("1" = 1, "2" = 1))

# Extract the coefficients of the hyperplane
coef <- t(svm_model$coefs) %*% svm_model$SV


# Update feature_map function to properly name the columns
feature_map <- function(X) {
	return(data.frame(X1 = X[,1], X2 = X[,2], Z = X[,1]^2 + X[,2]^2))
}

# Create a grid for the decision boundary
grid <- expand.grid(Var1 = seq(-3, 3, length.out = 100), Var2 = seq(-3, 3, length.out = 100))

# Map the grid to the feature space
Z_grid <- feature_map(grid)

names(grid)
head(grid)

head(Z_grid)
names(Z_grid)



# Compute Z value for the decision boundary
decision_boundary_surface <- -1*(coef[1]*Z_grid$X1 + coef[2]*Z_grid$X2 + svm_model$rho) / coef[3]

# Reshape decision boundary for plotting
Z_decision_boundary <- matrix(decision_boundary_surface, nrow = length(seq(-3, 3, length.out = 100)), ncol = length(seq(-3, 3, length.out = 100)))

# Plot the decision boundary in the transformed 3D space
plot_ly(data, x = ~Z1, y = ~Z2, z = ~Z3, color = ~y, colors = c("blue", "red"), type = "scatter3d", mode = "markers") %>%
	add_surface(x = Z_grid$X1, y = Z_grid$X2, z = Z_decision_boundary, showscale = FALSE)








# Map the decision boundary back to the original 2D space
decision_boundary_2d <- data.frame(X1 = grid$Var1, X2 = grid$Var2, Z = ifelse(Z_grid$decision_boundary > 0, 1, -1))


# Map the decision boundary back to the original 2D space
decision_boundary_2d <- data.frame(X1 = grid[,1], X2 = grid[,2], Z = ifelse(Z_grid$decision_boundary > 0, 1, -1))

# Plot the decision boundary in the original 2D space
ggplot(data, aes(x = X1, y = X2, color = as.factor(y))) +
	geom_point() +
	geom_contour(data = decision_boundary_2d, aes(x = X1, y = X2, z = Z), breaks = 0.5) +
	labs(x = "X1", y = "X2", title = "Original dataset with decision boundary") +
	theme_minimal()
