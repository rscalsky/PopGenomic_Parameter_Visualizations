library(ggplot2)
library(dplyr)
library(FSA)
library(segmented)
library(rpart)

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/R_inputs/")
getwd()

##### Data  #####

# PF3D7_0711200 #
API_df_rowwise <- read.csv("PF3D7_0711200_popgenomic_estimates_rowwise.csv")

# Removing country wide samples for countries with specific sites 0711200 ##

API_df_rowwise <- API_df_rowwise %>% ###remove this if generating graphs for all sites
  filter(!Geographic.Site %in% c("Cambodia", "Myanmar", "Senegal", "Vietnam"))


# PF3D7_0808100 #
API_df_rowwise <- read.csv("PF3D7_0808100_popgenomic_estimates_rowwise.csv")

#Removing country wide samples for countries with specific sites 0808100 #

API_df_rowwise <- API_df_rowwise %>% ###remove this if generating graphs for all sites
  filter(!Geographic.Site %in% c("Myanmar", "Senegal", "Vietnam"))



# PF3D7_1324300 #
API_df_rowwise <- read.csv("PF3D7_1324300_popgenomic_estimates_rowwise.csv")

#### Removing country wide samples for countries with specific sites 1324300 ##
API_df_rowwise <- API_df_rowwise %>% ###remove this if generating graphs for all sites
  filter(!Geographic.Site %in% c("Cambodia", "Mauritania", "Myanmar", "Senegal"))



# PF3D7_1335900 #
API_df_rowwise <- read.csv("PF3D7_1335900_popgenomic_estimates_rowwise.csv")

## Removing country wide samples for countries with specific sites 1335900 ##
API_df_rowwise <- API_df_rowwise %>% ###remove this if generating graphs for all sites
  filter(!Geographic.Site %in% c("Vietnam", "Myanmar", "Senegal"))



# 1361800 #
API_df_rowwise <- read.csv("PF3D7_1361800_popgenomic_estimates_rowwise.csv")

#### Removing country wide samples for countries with specific sites 1361800 ##
API_df_rowwise <- API_df_rowwise %>% ###remove this if generating graphs for all sites
  filter(!Geographic.Site %in% c("Cambodia", "Vietnam", "Myanmar", "Senegal"))

#### Analysis Start #####

#Keep only High API countries 
High_API_df_rowwise <- API_df_rowwise %>%
  filter(API_cat_adjusted == "High")

# Select Pi and API_regional variables 
clustering_df <- High_API_df_rowwise %>%
  select(Pi, API_regional)

# Select Pi and API_regional including moderate/low

clustering_all_df <- API_df_rowwise %>%
  select(Pi, API_regional)

#Scatter plot to visualize inflection point of Pi x API_regional 

ggplot(High_API_df_rowwise, aes(x = API_regional, y = Pi, color = Region)) +
  geom_point(size = 3) +
  labs(
    title = "Scatter Plot of Pi by API Regional",
    x = "API Regional",
    y = "Pi",
    color = "Region"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### K-means clustering on real data with testing for differences in groups #####

# k-means clustering on Pi values, cluster = 2
set.seed(123) # For reproducibility
clusters2 <- kmeans(High_API_df_rowwise$Pi, centers = 2)$cluster

# Add cluster = 2 assignments to the data frame
High_API_df_rowwise <- High_API_df_rowwise %>%
  mutate(Cluster2 = as.factor(clusters2))

# Visualize clusters = 2
ggplot(High_API_df_rowwise, aes(x = API_regional, y = Pi, color = Cluster2)) +
  geom_point(size = 3) +
  labs(
    title = "Scatter Plot of Pi by API Regional with Clusters",
    x = "API Regional",
    y = "Pi",
    color = "Cluster"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Statistical testing for differences (2 groups, k-cluster = 2)
# Wilcoxon rank-sum test
wilcox_test_result <- wilcox.test(Pi ~ Cluster2, data = High_API_df_rowwise)
print(wilcox_test_result)

# k-means clustering on Pi values, cluster = 3
set.seed(123) # For reproducibility
clusters3 <- kmeans(High_API_df_rowwise$API_regional, centers = 3)$cluster

# Add cluster = 3 assignments to the data frame
High_API_df_rowwise <- High_API_df_rowwise %>%
  mutate(Cluster3 = as.factor(clusters3))

# Visualize clusters = 3
ggplot(High_API_df_rowwise, aes(x = API_regional, y = Pi, color = Cluster3)) +
  geom_point(size = 3) +
  labs(
    title = "Scatter Plot of Pi by API Regional with Clusters",
    x = "API Regional",
    y = "Pi",
    color = "Cluster"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Perform Kruskal-Wallis test for 3 groups
kruskal_test_result <- kruskal.test(Pi ~ Cluster3, data = High_API_df_rowwise)
print(kruskal_test_result)

#post-hoc testing for Kruskal-wallis, Dunn test 
dunn_test_result <- dunnTest(Pi ~ Cluster3, data = High_API_df_rowwise, method = "bonferroni")
print(dunn_test_result)


# k-means clustering on two variables #

# Select Pi and API_regional variables for clustering
clustering_df <- High_API_df_rowwise %>%
  select(Pi, API_regional)

# Scale the data (standardization)
scaled_data <- scale(clustering_df)

# Set seed 
set.seed(123)

# Perform k-means clustering with 2 clusters
kmeans_result <- kmeans(scaled_data, centers = 3)

# Add cluster assignments to the original data frame
High_API_df_rowwise <- High_API_df_rowwise %>%
  mutate(Cluster = as.factor(kmeans_result$cluster))

# View clustering summary
print(kmeans_result)

# Scatter plot of clusters
ggplot(High_API_df_rowwise, aes(x = API_regional, y = Pi, color = Cluster)) +
  geom_point(size = 3) +
  labs(
    title = "K-Means Clustering of Pi and API Regional",
    x = "API Regional",
    y = "Pi",
    color = "Cluster"
  ) +
  theme_minimal()

# Determine the optimal number of clusters using the Elbow Method
wss <- sapply(1:10, function(k) {
  kmeans(scale(clustering_df), centers = k, nstart = 20)$tot.withinss
})

# Plot WSS against the number of clusters
plot(1:10, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters (k)", ylab = "Total Within-Cluster Sum of Squares",
     main = "Elbow Method for Optimal k")


#### Derivative reate of change approach real data ####

#////////////////////////////////////////

# Derivative (Rate of Change) Approach 

# Calculate the rate of change (derivative) of Pi with respect to API_regional
clustering_df <- clustering_df %>%
  arrange(API_regional) %>%
  mutate(rate_of_change = c(NA, diff(Pi) / diff(API_regional)))

# Plot the rate of change
ggplot(clustering_df, aes(x = API_regional, y = rate_of_change)) +
  geom_point() +
  labs(title = "Rate of Change of Pi vs API_regional", x = "API_regional", y = "Rate of Change")

# Find the point with the maximum rate of change
max_rate_point <- clustering_df[which.max(clustering_df$rate_of_change), ]
max_rate_point

# Top 3 rate of change 
# Sort the clustering_df by rate_of_change in descending order and get the top 3 points
top_3_rate_points <- clustering_df %>%
  arrange(desc(rate_of_change)) %>%
  head(3)

# View the top 3 points with the highest rate of change
top_3_rate_points

#### Simulated Analyses #####

#////////////////////////////////////////

# Fit smoothed model # 

# Fit a smooth, non-linear line using loess (locally weighted regression)
loess_model <- loess(Pi ~ API_regional, data = clustering_df)

# Plot the original data and the fitted smooth line
ggplot(clustering_df, aes(x = API_regional, y = Pi)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, col = "blue") +
  labs(title = "Smooth Line of Best Fit (Loess)", x = "API_regional", y = "Pi")



# Generate simulated Pi values for every 10 API 

# Create a sequence of API_regional values from the min to max of the data
api_sequence <- seq(min(clustering_df$API_regional), max(clustering_df$API_regional), by = 10)

# Predict Pi values for the API_regional sequence
predicted_Pi <- predict(loess_model, newdata = data.frame(API_regional = api_sequence))

# Combine API_regional and predicted Pi values into a new dataframe
predicted_df <- data.frame(API_regional = api_sequence, Pi = predicted_Pi)

# View the first few predicted values
head(predicted_df)



# Calculate the rate of change (derivative) for the predicted data
predicted_df <- predicted_df %>%
  mutate(rate_of_change = c(NA, diff(Pi) / diff(API_regional)))

# Plot the rate of change
ggplot(predicted_df, aes(x = API_regional, y = rate_of_change)) +
  geom_point() +
  labs(title = "Rate of Change of Pi (Predicted)", x = "API_regional", y = "Rate of Change")

# Sort the predicted data by rate_of_change in descending order and get the top 3 points
top_3_predicted_rate_points <- predicted_df %>%
  arrange(desc(rate_of_change)) %>%
  head(3)

# View the top 3 points with the highest rate of change
top_3_predicted_rate_points




# K-means clustering on rates of change in simulated rate of change

# Select the rate_of_change variable for clustering
clustering_df_predicted <- predicted_df %>%
  select(rate_of_change, API_regional)

# Remove rows with NA values in the rate_of_change column
clustering_df_predicted_clean <- clustering_df_predicted %>%
  filter(!is.na(rate_of_change))

# Perform k-means clustering (e.g., 3 clusters)
set.seed(123)  # Set a seed for reproducibility
kmeans_result <- kmeans(clustering_df_predicted_clean, centers = 6)


# Remove rows with NA values in the rate_of_change column
predicted_df_clean <- predicted_df %>%
  filter(!is.na(rate_of_change))

# Add the cluster assignment to the predicted_df
predicted_df_clean$cluster <- kmeans_result$cluster

# View the k-means result
head(predicted_df)

# Plot the rate of change with k-means clusters
ggplot(predicted_df_clean, aes(x = API_regional, y = rate_of_change, color = factor(cluster))) +
  geom_point() +
  labs(title = "K-means Clustering of Rate of Change", x = "API_regional", y = "Rate of Change", color = "Cluster") +
  scale_color_manual(values = c("red", "blue", "green", "orange", "black", "purple"))  # Custom colors for clusters
  




# K-means clustering on rates of change in simulated Pi

# Select the rate_of_change variable for clustering
clustering_df_predicted_Pi <- predicted_df %>%
  select(Pi, API_regional)

# Perform k-means clustering (e.g., 3 clusters)
set.seed(123)  # Set a seed for reproducibility
kmeans_result <- kmeans(clustering_df_predicted_Pi, centers = 4)

# Add the cluster assignment to the predicted_df
predicted_df$cluster <- kmeans_result$cluster

# View the k-means result
head(predicted_df)

# Plot the rate of change with k-means clusters
ggplot(predicted_df, aes(x = API_regional, y = Pi, color = factor(cluster))) +
  geom_point() +
  labs(title = "K-means Clustering of Simulated Pi", x = "API_regional", y = "Pi", color = "Cluster") +
  scale_color_manual(values = c("red", "blue", "green", "orange"))  # Custom colors for clusters


#### Alternative methods #####

#////////////////////////////////////////

# Fit a decision tree
tree_model <- rpart(Pi ~ API_regional, data = predicted_df, control = rpart.control(minisplit = 2, cp = 0.001))

# Visualize the tree
plot(tree_model)
text(tree_model, use.n = TRUE)

# Print the tree structure
print(tree_model)


#////////////////////////////////////////


# Segmented regression approach #

# Fit a linear model
lm_model <- lm(Pi ~ API_regional, data = High_API_df_rowwise)

# Apply segmented regression to find the breakpoint
seg_model <- segmented(lm_model, seg.Z = ~ API_regional)

# View the results
summary(seg_model)

#extract fitted values 
fitted_values <- fitted(seg_model)

#df for plotting 
segmented_plot_data <- data.frame(API_regional = High_API_df_rowwise$API_regional,
                                  Pi = High_API_df_rowwise$Pi,
                                  Fitted_Pi = fitted_values)

# Plot the segmented regression
plot(High_API_df_rowwise$API_regional, High_API_df_rowwise$Pi, pch = 19, col = "blue",
     xlab = "API Regional", ylab = "Pi")
lines(seg_model, col = "red")


ggplot(segmented_plot_data, aes(x = API_regional, y = Pi)) +
  geom_point() +  # Scatter plot of the original data
  geom_line(aes(x = API_regional, y = Fitted_Pi), color = "blue") +  # Segmented regression line
  labs(title = "Segmented Regression of Pi vs API_regional", 
       x = "API_regional", 
       y = "Pi")

#extract breakpoint APIs

breakpoints <- seg_model$psi

print(breakpoints)