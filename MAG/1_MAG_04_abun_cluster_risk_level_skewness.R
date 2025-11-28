MAG_abun <- read.csv("/path/to/data/source_data_1.csv", row.names = 1)
data_s <- read.csv("/path/to/data/source_data_2.csv", row.names = 1)

# Extract MAG names based on cluster (CR) type
MAG_names_A <- rownames(data_s %>% filter(CR == "Class A"))
MAG_names_B <- rownames(data_s %>% filter(CR == "Class B"))
MAG_names_C <- rownames(data_s %>% filter(CR == "Class C"))

# Subset MAG abundance for Class A and calculate row sums
MAG_abun_A <- MAG_abun[MAG_names_A, ] %>%
  t() %>%
  data.frame() %>%
  mutate(sum = rowSums(.))

# Load necessary packages
library(stats) # for kmeans function
library(Rtsne) # for t-SNE algorithm
library(factoextra) # for k-means evaluation (optimal k)
library(cluster) # for silhouette score calculation
set.seed(13)

# Calculate skewness for sum values and log-transformed sum values
sum_values <- MAG_abun_A$sum
library(moments)
skewness(sum_values) # skewness for original sum values
skewness(log(sum_values + 0.00000001)) # skewness for log-transformed sum values
skewness(log(sum_values + 1)) # skewness with +1 log transformation
y_values <- log(sum_values + 1) # log transformation (log(x+1))

# Perform Kolmogorov-Smirnov test for normality
ks.test(y_values, "pnorm", mean = mean(y_values), sd = sd(y_values))

# Create histogram with normal distribution curve for log-transformed data
p <- ggplot(data.frame(y = y_values), aes(x = y)) +
  geom_histogram(aes(y = ..density..),
    bins = 30, # Adjust the number of bins
    color = "black",
    fill = "lightblue"
  ) +
  stat_function(
    fun = dnorm,
    args = list(mean = mean(y_values), sd = sd(y_values)),
    color = "red",
    size = 1.2
  ) +
  geom_vline(
    xintercept = mean(y_values),
    linetype = "dashed",
    color = "blue"
  ) +
  labs(
    x = "ln(ARG abundance + 1)", # Label for x-axis
    y = "Distribution density" # Label for y-axis
  ) +
  theme_base() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "darkgray"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
ggsave(p, filename = "/path/to/figure/distribution_density_logxplus1.pdf")

# Perform KS test on original sum values for normality
ks.test(sum_values, "pnorm", mean = mean(sum_values), sd = sd(sum_values))

# Create histogram with normal distribution curve for original sum values
p2 <- ggplot(data.frame(y = sum_values), aes(x = y)) +
  geom_histogram(aes(y = ..density..),
    bins = 30, # Adjust the number of bins
    color = "black",
    fill = "lightblue"
  ) +
  stat_function(
    fun = dnorm,
    args = list(mean = mean(sum_values), sd = sd(sum_values)),
    color = "red",
    size = 1.2
  ) +
  geom_vline(
    xintercept = mean(sum_values),
    linetype = "dashed",
    color = "blue"
  ) +
  labs(
    x = "ARG abundance", # Label for x-axis
    y = "Distribution density" # Label for y-axis
  ) +
  theme_base() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "darkgray"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
ggsave(p2, filename = "/path/to/figure/distribution_density_x.pdf")

# Create Q-Q plot for log-transformed data
pdf("/path/to/figure/res_figure_3.pdf", width = 6, height = 6)
qqnorm(y_values)
qqline(y_values, col = "red")
dev.off()

# Create Q-Q plot for original sum values
pdf("/path/to/figure/res_figure_4.pdf", width = 6, height = 6)
qqnorm(sum_values)
qqline(sum_values, col = "red")
dev.off()

# Calculate minimum sample size based on standard deviation and confidence level
sigma <- sd(y_values)
confidence_level <- 0.99 # 99% confidence level
delta <- 0.05 # Allowable margin of error
z_value <- qnorm(1 - (1 - confidence_level) / 2) # Z-value for two-tailed normal distribution
n_min <- ((z_value * sigma) / delta)^2 # Minimum sample size
cat("Minimum sample size (original):", round(n_min, 0), "\n")

# Adjust for dropout rate (20%)
dropout_rate <- 0.2
n_adjusted <- n_min / (1 - dropout_rate)
cat("Minimum sample size after dropout adjustment:", round(n_adjusted, 0), "\n")

# Step 1: Find the optimal k-value using silhouette scores and other metrics
k_values <- 4:16
library(clusterSim)
silhouette_scores <- numeric(length(k_values))
davies_bouldin_scores <- numeric(length(k_values))
calinski_harabasz_scores <- numeric(length(k_values))

# Loop through different k-values to calculate evaluation scores
for (i in seq_along(k_values)) {
  k <- k_values[i]

  # Perform k-means clustering
  kmeans_result <- kmeans(y_values, centers = k, nstart = 25)

  # Calculate Silhouette Score (higher is better)
  silhouette_avg <- mean(silhouette(kmeans_result$cluster, dist(as.matrix(y_values)))[, 3])
  silhouette_scores[i] <- silhouette_avg

  # Calculate Davies-Bouldin Score (lower is better)
  davies_bouldin_scores[i] <- index.DB(as.matrix(y_values), cl = kmeans_result$cluster)$DB

  # Calculate Calinski-Harabasz Score (higher is better)
  cluster_stats <- cluster.stats(d = dist(as.matrix(y_values)), clustering = kmeans_result$cluster)
  calinski_harabasz_scores[i] <- cluster_stats$ch
}

# Save results of the evaluation
results <- data.frame(
  k = k_values,
  Silhouette_Score = silhouette_scores,
  Davies_Bouldin_Score = davies_bouldin_scores,
  Calinski_Harabasz_Score = calinski_harabasz_scores
)

# Step 2: Plot silhouette, Davies-Bouldin, and Calinski-Harabasz scores
library(patchwork)
library(ggthemes)
library(ggplot2)

font_size <- 10
axis_font_size <- 8
x_value <- 10 # Adjust based on the best k-value found

# Plot silhouette score
p1 <- ggplot(results, aes(x = k, y = Silhouette_Score)) +
  geom_line(color = "orange") +
  geom_point(color = "orange3") +
  geom_vline(xintercept = x_value, linetype = "dashed", color = "black") +
  labs(title = "Silhouette Score", x = "Number of Clusters (K)", y = "Silhouette Score") +
  theme_base() +
  theme(
    plot.title = element_text(size = font_size, hjust = 0.5),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_font_size)
  )

# Plot Davies-Bouldin score
p2 <- ggplot(results, aes(x = k, y = Davies_Bouldin_Score)) +
  geom_line(color = "orange") +
  geom_point(color = "orange3") +
  geom_vline(xintercept = x_value, linetype = "dashed", color = "black") +
  labs(title = "Davies-Bouldin Score", x = "Number of Clusters (K)", y = "Davies-Bouldin Score") +
  theme_base() +
  theme(
    plot.title = element_text(size = font_size, hjust = 0.5),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_font_size)
  )

# Plot Calinski-Harabasz score
p3 <- ggplot(results, aes(x = k, y = Calinski_Harabasz_Score)) +
  geom_line(color = "orange") +
  geom_point(color = "orange3") +
  geom_vline(xintercept = x_value, linetype = "dashed", color = "black") +
  labs(title = "Calinski-Harabasz Score", x = "Number of Clusters (K)", y = "Calinski - Harabasz Score") +
  theme_base() +
  theme(
    plot.title = element_text(size = font_size, hjust = 0.5),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_font_size)
  )

# Use patchwork to combine the three plots into one layout
combined_plot <- p1 + p2 + p3 + plot_layout(nrow = 1)
combined_plot

# Save the combined plot to a PDF file
ggsave(combined_plot,
  filename = "/path/to/figure/best_kmeans_value_three_scores.pdf",
  family = "ArialMT", height = 4, width = 12
)

# Select the best k value based on silhouette score
best_k_silhouette <- 10 # Assume the best k value based on silhouette score is 10

# Step 2: Perform k-means clustering with the best k value
set.seed(13)
kmeans_result <- kmeans(y_values, centers = best_k_silhouette, nstart = 100, iter.max = 100)

# Add the clustering results to the MAG abundance dataframe
MAG_abun_A$cluster <- as.factor(kmeans_result$cluster)

# Reorder clusters by the mean of sum values
cluster_means <- tapply(y_values, MAG_abun_A$cluster, mean)
sorted_clusters <- order(cluster_means)

# Create a mapping to assign new cluster labels
cluster_mapping <- setNames(paste0("V", 1:best_k_silhouette), sorted_clusters)
MAG_abun_A$cluster <- cluster_mapping[as.character(MAG_abun_A$cluster)]

# Output the frequency of each cluster
MAG_abun_A$cluster %>% table()

# Step 3: Use t-SNE for 2D visualization of the clustering results
# t-SNE input requires a matrix format
input_data <- as.matrix(y_values)
tsne_result <- Rtsne(input_data, dims = 2, check_duplicates = FALSE)

# Create a new dataframe for the t-SNE results
tsne_df <- data.frame(
  x = tsne_result$Y[, 1],
  y = tsne_result$Y[, 2],
  cluster = MAG_abun_A$cluster
)

# Step 4: Plot the t-SNE results in 2D
p_tsne <- ggplot(tsne_df, aes(x = x, y = y, color = cluster)) +
  geom_point() +
  labs(
    title = "t-SNE visualization of k-means clustering",
    x = "t-SNE dimension 1",
    y = "t-SNE dimension 2"
  ) +
  theme_minimal()

# Save the t-SNE plot to a PDF file
ggsave(p_tsne, filename = "/path/to/figure/kmeans_tsne_graph.pdf", width = 6, height = 6)

# Step 5: Create a function to predict the cluster of new data points
predict_kmeans <- function(new_data, centers) {
  # Compute the distance from each new data point to the cluster centers
  distances <- apply(centers, 1, function(center) {
    sqrt(rowSums((new_data - center)^2)) # Euclidean distance
  })
  # Return the index of the closest cluster center
  max.col(-distances) # Find the column with the smallest distance
}

# If the original data is a single-column matrix/data frame, convert new data to the same structure
new_data_matrix <- matrix(y_values, ncol = 1)

# Predict the cluster labels for the new data
predicted_clusters <- predict_kmeans(new_data_matrix, kmeans_result$centers)

# Print predicted clusters along with the original cluster labels
data.frame(pred = predicted_clusters, obs = kmeans_result$cluster)

# Save the final data with cluster information to a CSV file
MAG_abun_A %>% write.csv("/path/to/figure/ML/models_y_1_original.csv")
