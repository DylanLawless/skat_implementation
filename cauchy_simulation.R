# cauchy transformation of P-values ----
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare the data
set.seed(123)
n <- 100  # Number of P-values
p_values <- runif(n, 0, 1)
cauchy_transformed <- tan(pi * (0.5 - p_values))
data <- data.frame(Index = 1:n, P_Values = p_values, Transformed = cauchy_transformed)

# Convert from wide to long format
data_long <- data %>%
  pivot_longer(cols = c("P_Values", "Transformed"), 
               names_to = "Type", 
               values_to = "Value")

# View the first few rows of the long format data
head(data_long)

p_dotplots_log <- 
  data_long |>
  filter(Type == "P_Values") |>
  ggplot(aes(x = Index, y = -log10(Value), color = Type)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Original P-values",
       x = "Index",
       y = "-log10(Values)") +
  scale_color_manual(values = c("P_Values" = "blue", "Transformed" = "red"))

print(p_dotplots_log)
ggsave("./images/p_cauchy_dotplots_log.pdf", plot = p_dotplots_log)

# Histograms
p_histograms <- ggplot(data_long, aes(x = Value, fill = Type)) +
  geom_histogram(bins = 20, alpha = 0.5) +
  facet_wrap(~ Type, scales = "free_x") +
  theme_minimal() +
  labs(title = "Original and Transformed P-values",
       x = "Values",
       y = "Frequency") +
  scale_fill_manual(values = c("P_Values" = "blue", "Transformed" = "red"))

print(p_histograms)
ggsave("./images/p_cauchy_histograms.pdf", plot = p_histograms)

# Dotplots
p_dotplots <- ggplot(data_long, aes(x = Index, y = Value, color = Type)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Original and Transformed P-values",
       x = "Index",
       y = "Values") +
  scale_color_manual(values = c("P_Values" = "blue", "Transformed" = "red"))

print(p_dotplots)
ggsave("./images/p_cauchy_dotplots.pdf", plot = p_dotplots)

# Scatter plot of P-values vs. Transformed Values
p_scatter <- ggplot(data, aes(x = P_Values, y = Transformed)) +
  geom_point(alpha = 0.7, color = "blue") +
  geom_line() +
  theme_minimal() +
  labs(title = "Original P-values vs. Transformed Values",
       x = "Original P-values",
       y = "Transformed Values")

# Display the plot
print(p_scatter)
ggsave("./images/p_cauchy_scatter.pdf", plot = p_scatter)

# Tan function ----

# Load necessary library
library(ggplot2)

# Create a sequence of values for x (choosing a range around 0 to highlight interesting parts)
x_values <- seq(-2 * pi, 2 * pi, length.out = 1000)

# Apply the tangent function
y_values <- tan(x_values)

# Create a dataframe for plotting
tan_data <- data.frame(X = x_values, Y = y_values)

# Plot the tangent function
tan_plot <- ggplot(tan_data, aes(x = X, y = Y)) +
  geom_line(color = "blue") +
  theme_minimal() +
  labs(title = "Tangent Function (tan(x))",
       x = "x",
       y = "tan(x)") +
  ylim(-10, 10) # Limiting the y-axis to avoid extreme values

# Display the plot
print(tan_plot)
ggsave("./images/p_cauchy_tan_plot.pdf", plot = tan_plot)

# Pi term ----

# Load necessary library
library(ggplot2)

# Create a sequence of P-values from 0 to 1
p_values <- seq(0, 1, length.out = 1000)

# Apply the transformation with and without the pi term
transformed_with_pi <- tan(pi * (0.5 - p_values))
transformed_without_pi <- tan(0.5 - p_values)

# Create a dataframe for plotting
transform_data <- data.frame(P_Values = p_values, 
                             Transformed_with_Pi = transformed_with_pi, 
                             Transformed_without_Pi = transformed_without_pi)

# Convert data from wide to long format
transform_data_long <- tidyr::pivot_longer(transform_data, 
                                           cols = c("Transformed_with_Pi", "Transformed_without_Pi"), 
                                           names_to = "Transformation", 
                                           values_to = "Value")

# Plot the transformations
pi_effect_plot <- ggplot(transform_data_long, aes(x = P_Values, y = Value, color = Transformation)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Effect of Pi in Cauchy Distribution Transformation",
       x = "Original P-values",
       y = "Transformed Values") +
  scale_color_manual(values = c("Transformed_with_Pi" = "blue", 
                                "Transformed_without_Pi" = "red"))

# Display the plot
print(pi_effect_plot)
ggsave("./images/p_cauchy_pi_effect_plot.pdf", plot = pi_effect_plot)

# Pi visualised ----

# Load necessary library
library(ggplot2)

# Define a circle
t <- seq(0, 2*pi, length.out = 100)
x <- cos(t)
y <- sin(t)
circle_data <- data.frame(x = x, y = y)

# Plot the circle with diameter, circumference, and area annotations
combined_circle_plot <- ggplot(circle_data, aes(x, y)) +
  geom_polygon(fill = "lightblue", color = "blue", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0), arrow = arrow(type = "closed"))+
  annotate("text", x = 0.5, y = 0.2, label = "Diameter = 1")+
  annotate("text", x = 0, y = -1.2, label = "Circumference = π * Diameter") +
  annotate("text", x = 0, y = 1.2, label = "Area = π * r^2 (r = 0.5)")+
  coord_fixed() +
  theme_void() +
  labs(title = "Properties of a Circle Related to π")

# Display the plot
print(combined_circle_plot)
ggsave("./images/p_cauchy_combined_circle_plot.pdf", plot = combined_circle_plot)

# Why pi ? ---

# Load necessary library
library(ggplot2)

# Create a sequence of P-values
p_values <- seq(0, 1, length.out = 1000)

# Apply transformations using pi and an arbitrary number (e.g., 3)
transformed_with_pi <- tan(pi * (0.5 - p_values))
transformed_with_3 <- tan(3 * (0.5 - p_values))

# Create a dataframe for plotting
comparison_data <- data.frame(
  P_Values = p_values,
  Transformed_with_Pi = transformed_with_pi,
  Transformed_with_3 = transformed_with_3
)

# Convert data from wide to long format
comparison_data_long <- tidyr::pivot_longer(comparison_data, 
                                            cols = c("Transformed_with_Pi", "Transformed_with_3"), 
                                            names_to = "Transformation", 
                                            values_to = "Value")

# Plot the transformations
comparison_plot <- ggplot(comparison_data_long, aes(x = P_Values, y = Value, color = Transformation)) +
  geom_line() +
  geom_line(data = subset(comparison_data_long, Transformation == "Transformed_with_3"),
            linetype = "dotted", size = 1) +
  theme_minimal() +
  labs(title = "Comparison of Transformations\nUsing Pi vs. Arbitrary Number",
       x = "Original P-values",
       y = "Transformed Values") +
  scale_color_manual(values = c("Transformed_with_Pi" = "blue", 
                                "Transformed_with_3" = "red"))

# Display the plot
comparison_plot
ggsave("./images/p_cauchy_cauchy_comparison_plot.pdf", plot = comparison_plot)







