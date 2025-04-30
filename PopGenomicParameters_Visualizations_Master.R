library(ggpattern) 
library(ggplot2)
library(dplyr)

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/R_inputs/")
getwd()

PF3D7_0711200_rowwise <- read.csv("PF3D7_0711200_popgenomic_estimates_rowwise.csv")

##### Nucleotide Diversity #####

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat <- factor(PF3D7_0711200_rowwise$API_cat, levels = c("High", "Moderate", "Low"))

# Arrange the data by API_cat and descending Pi
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(API_cat) %>%
  arrange(desc(Pi), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot with error bars
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = Pi, fill = API_cat)) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) + 
  geom_errorbar(aes(ymin = Pi - PiSD, ymax = Pi + PiSD), 
                width = 0.2, 
                position = position_dodge(0.7)) + 
  labs(x = "Geographic Site", y = "Pi", fill = "API Category") +
  ggtitle("Global Nucleotide Diversity of PF3D7_0711200") +  # Add title here
    theme_minimal() +
  scale_fill_manual(values = c("High" = "darkred", "Moderate" = "darkorange", "Low" = "darkgreen")) +
  theme(
    plot.title = element_text(hjust = 0.5), #center the title 
    axis.text.x = element_text(angle = 45, hjust = 1)) # Optional: tilt x-axis labels for readability

##### Nucleotide Diversity With shading for region ######

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat <- factor(PF3D7_0711200_rowwise$API_cat, levels = c("High", "Moderate", "Low"))

# Convert Region to a factor with the desired order
PF3D7_0711200_rowwise$Region <- factor(PF3D7_0711200_rowwise$Region, 
                                       levels = c("Central Africa", "East Africa", "West Africa", "South East Asia"))

# Arrange the data by Region, API_cat, and descending Pi
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(Region, API_cat) %>%
  arrange(desc(Pi), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot with error bars and distinguish regions with patterns
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = Pi_pct, fill = API_cat, pattern = Region)) + 
  geom_bar_pattern(stat = "identity", position = position_dodge(), width = 0.7, 
                   pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "white") + 
  geom_errorbar(aes(ymin = Pi_pct - PiSD_pct, ymax = Pi_pct + PiSD_pct), 
                width = 0.2, 
                position = position_dodge(0.7)) + 
  labs(x = "Geographic Site", y = "Pi (%)", fill = "API Category", pattern = "Region") +
  ggtitle("Global Nucleotide Diversity of Conserved Pf Protein, Unknown Function (PF3D7_0711200)") +  # Add title here
  scale_y_continuous(breaks = seq(0, 0.65, by = .05), limits= c(0, 0.65)) +
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +
  scale_pattern_manual(values = c("Central Africa" = "stripe", 
                                  "East Africa" = "crosshatch", 
                                  "West Africa" = "wave",  # You can set horizontal stripes here
                                  "South East Asia" = "circle")) + # or any other pattern
  theme(
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size 
    plot.title = element_text(hjust = 0.5),  # Center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Optional: tilt x-axis labels for readability
    axis.text.y = element_text(size =12)
  )


##### Haplotype Diversity ####

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat <- factor(PF3D7_0711200_rowwise$API_cat, levels = c("High", "Moderate", "Low"))

# Arrange the data by API_cat and descending Pi
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(API_cat) %>%
  arrange(desc(HaplotypeDiversity), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot with error bars
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = HaplotypeDiversity, fill = API_cat, pattern = Region)) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) + 
  geom_errorbar(aes(ymin = HaplotypeDiversity - HaplotypeSD, ymax = HaplotypeDiversity + HaplotypeSD), 
                width = 0.2, 
                position = position_dodge(0.7)) + 
  labs(x = "Geographic Site", y = "Haplotype Diversity", fill = "API Category", pattern = "Region") +
  ggtitle("Global Haplotype Diversity of PF3D7_0711200") +  # Add title here
  theme_minimal() +
  scale_fill_manual(values = c("High" = "darkred", "Moderate" = "darkorange", "Low" = "darkgreen")) +
  theme(
    plot.title = element_text(hjust = 0.5), #center the title 
    axis.text.x = element_text(angle = 45, hjust = 1)) # Optional: tilt x-axis labels for readability

##### Haplotype Diversity with Shading for Region #####

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat <- factor(PF3D7_0711200_rowwise$API_cat, levels = c("High", "Moderate", "Low"))

# Convert Region to a factor with the desired order
PF3D7_0711200_rowwise$Region <- factor(PF3D7_0711200_rowwise$Region, 
                                       levels = c("Central Africa", "East Africa", "West Africa", "South East Asia"))

# Arrange the data by API_cat and descending Pi
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(Region, API_cat) %>%
  arrange(desc(HaplotypeDiversity), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot with error bars
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = HaplotypeDiversity, fill = API_cat, pattern = Region)) + 
  geom_bar_pattern(stat = "identity", position = position_dodge(), width = 0.7,
           pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "white") + 
  geom_errorbar(aes(ymin = HaplotypeDiversity - HaplotypeSD, ymax = HaplotypeDiversity + HaplotypeSD), 
                width = 0.2, 
                position = position_dodge(0.7)) + 
  labs(x = "Geographic Site", y = "Haplotype Diversity", fill = "API Category", pattern = "Region") +
  ggtitle("Global Haplotype Diversity of Conserved Pf Protein, Unknown Function (PF3D7_0711200)") +  # Add title here
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +
  scale_pattern_manual(values = c("Central Africa" = "stripe", 
                                  "East Africa" = "crosshatch", 
                                  "West Africa" = "wave",  # You can set horizontal stripes here
                                  "South East Asia" = "circle")) + # or any other pattern
  theme(
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size   
    plot.title = element_text(hjust = 0.5), #center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size =12)) # Optional: tilt x-axis labels for readability

##### Number of Haplotypes ####

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat <- factor(PF3D7_0711200_rowwise$API_cat, levels = c("High", "Moderate", "Low"))

# Arrange the data by API_cat and descending Pi
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(API_cat) %>%
  arrange(desc(NumberHaplotypes), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot with error bars
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = NumberHaplotypes, fill = API_cat)) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) + 
  labs(x = "Geographic Site", y = "Number of Haplotypes", fill = "API Category") +
  ggtitle("Global Haplotype Counts of PF3D7_0711200") +  # Add title here
  theme_minimal() +
  scale_fill_manual(values = c("High" = "darkred", "Moderate" = "darkorange", "Low" = "darkgreen")) +
  theme(
    plot.title = element_text(hjust = 0.5), #center the title 
    axis.text.x = element_text(angle = 45, hjust = 1)) # Optional: tilt x-axis labels for readability

##### Number of Haplotypes with Region shading #####

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat <- factor(PF3D7_0711200_rowwise$API_cat, levels = c("High", "Moderate", "Low"))

# Convert Region to a factor with the desired order
PF3D7_0711200_rowwise$Region <- factor(PF3D7_0711200_rowwise$Region, 
                                       levels = c("Central Africa", "East Africa", "West Africa", "South East Asia"))

# Arrange the data by API_cat and descending Pi
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(Region, API_cat) %>%
  arrange(desc(NumberHaplotypes), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot with error bars
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = NumberHaplotypes, fill = API_cat, pattern = Region)) + 
  geom_bar_pattern(stat = "identity", position = position_dodge(), width = 0.7,
                   pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "white") + 
  labs(x = "Geographic Site", y = "Number of Haplotypes", fill = "API Category", pattern = "Region") +
  ggtitle("Global Haplotype Counts of Conserved Pf Protein, Unknown Function (PF3D7_0711200)") +  # Add title here
  scale_y_continuous(breaks = seq(0, 150, by = 10), limits= c(0, 150)) +
  scale_pattern_manual(values = c("Central Africa" = "stripe", 
                                  "East Africa" = "crosshatch", 
                                  "West Africa" = "wave",  # You can set horizontal stripes here
                                  "South East Asia" = "circle")) + # or any other pattern
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +
  theme(
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size
    plot.title = element_text(hjust = 0.5), #center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)) # Optional: tilt x-axis labels for readability

##### Number of Segregating Sites #####

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat <- factor(PF3D7_0711200_rowwise$API_cat, levels = c("High", "Moderate", "Low"))

# Arrange the data by API_cat and descending Pi
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(API_cat) %>%
  arrange(desc(NumberSegregatingSites), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot with error bars
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = NumberSegregatingSites, fill = API_cat)) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) + 
  labs(x = "Geographic Site", y = "Numberof Segregating Sites", fill = "API Category") +
  ggtitle("Global Segregating Site Counts of PF3D7_0711200") +  # Add title here
  theme_minimal() +
  scale_fill_manual(values = c("High" = "darkred", "Moderate" = "darkorange", "Low" = "darkgreen")) +
  theme(
    plot.title = element_text(hjust = 0.5), #center the title 
    axis.text.x = element_text(angle = 45, hjust = 1)) # Optional: tilt x-axis labels for readability

##### Number of Segregating Sites with Region Shading #####

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat <- factor(PF3D7_0711200_rowwise$API_cat, levels = c("High", "Moderate", "Low"))

# Convert Region to a factor with the desired order
PF3D7_0711200_rowwise$Region <- factor(PF3D7_0711200_rowwise$Region, 
                                       levels = c("Central Africa", "East Africa", "West Africa", "South East Asia"))

# Arrange the data by API_cat and descending Pi
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(Region, API_cat) %>%
  arrange(desc(NumberSegregatingSites), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot with error bars
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = NumberSegregatingSites, fill = API_cat, pattern = Region)) + 
  geom_bar_pattern(stat = "identity", position = position_dodge(), width = 0.7,
           pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "white") + 
  labs(x = "Geographic Site", y = "Numberof Segregating Sites", fill = "API Category", "Region") +
  ggtitle("Global Segregating Site Counts of PF3D7_0711200") +  # Add title here
  scale_pattern_manual(values = c("Central Africa" = "stripe", 
                                  "East Africa" = "crosshatch", 
                                  "West Africa" = "wave",  # You can set horizontal stripes here
                                  "South East Asia" = "circle")) + # or any other pattern
  theme_minimal() +
  scale_fill_manual(values = c("High" = "darkred", "Moderate" = "darkorange", "Low" = "darkgreen")) +
  theme(
    plot.title = element_text(hjust = 0.5), #center the title 
    axis.text.x = element_text(angle = 45, hjust = 1)) # Optional: tilt x-axis labels for readability

##### Tajima's D #####

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat <- factor(PF3D7_0711200_rowwise$API_cat, levels = c("High", "Moderate", "Low"))

# Arrange the data by API_cat and Tajima'sD (descending order for consistency)
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(API_cat) %>%
  arrange(desc(Tajima.sD), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot for Tajima'sD
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = Tajima.sD, fill = API_cat)) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) + 
  geom_hline(yintercept = 2, linetype = "dotted", color = "red") +  # Add dotted red line at y = 2.5
  geom_hline(yintercept = -2, linetype = "dotted", color = "red") + # Add dotted red line at y = -2.5
    labs(x = "Geographic Site", y = "Tajima's D", fill = "API Category") +
  ggtitle("Global Tajima's D of PF3D7_0711200") +  # Add title here
  ylim(-2.5, 2.5) +  # Set y-axis limits
  theme_minimal() +
  scale_fill_manual(values = c("High" = "darkred", "Moderate" = "darkorange", "Low" = "darkgreen")) +
  theme(
    plot.title = element_text(hjust = 0.5),    # Center the title
    axis.text.x = element_text(angle = 45, hjust = 1) # Optional: tilt x-axis labels for readability
  )

##### Tajima's D with Region Shading #####

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat <- factor(PF3D7_0711200_rowwise$API_cat, levels = c("High", "Moderate", "Low"))

# Convert Region to a factor with the desired order
PF3D7_0711200_rowwise$Region <- factor(PF3D7_0711200_rowwise$Region, 
                                       levels = c("Central Africa", "East Africa", "West Africa", "South East Asia"))

# Arrange the data by API_cat and Tajima'sD (descending order for consistency)
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(Region, API_cat) %>%
  arrange(desc(Tajima.sD), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot for Tajima'sD
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = Tajima.sD, fill = API_cat, pattern = Region)) + 
  geom_bar_pattern(stat = "identity", position = position_dodge(), width = 0.7,
                   pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "white") + 
  geom_hline(yintercept = 2, linetype = "dotted", color = "red", size = 1.4) +  # Add dotted red line at y = 2.5
  geom_hline(yintercept = -2, linetype = "dotted", color = "red", size = 1.4) + # Add dotted red line at y = -2.5
  labs(x = "Geographic Site", y = "Tajima's D", fill = "API Category", pattern = "Region") +
  ggtitle("Global Tajima's D of Conserved Pf Protein, Unknown Function (PF3D7_0711200)") +  # Add title here
  scale_pattern_manual(values = c("Central Africa" = "stripe", 
                                  "East Africa" = "crosshatch", 
                                  "West Africa" = "wave",  # You can set horizontal stripes here
                                  "South East Asia" = "circle")) + # or any other pattern
  scale_y_continuous(breaks = seq(-3, 3, by = 0.5), limits= c(-3,3)) +
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +
  theme(
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size   
    plot.title = element_text(hjust = 0.5), #center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size =12) # Optional: tilt x-axis labels for readability
  )

##### Regressions #####

# Perform the linear regression
model <- lm(Pi_pct ~ API, data = PF3D7_0711200_rowwise)

# Display the summary of the regression model
summary(model)

coefficients <- summary(model)$coefficients

slope <- coefficients["API", "Estimate"]
intercept <- coefficients["(Intercept)", "Estimate"]

formula <- paste("Pi = ", slope, " * API + ", intercept)
print(formula)

##### Dot plot of Pi x API with line of best fit #####
library(broom)  # For tidy regression output

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat <- factor(PF3D7_0711200_rowwise$API_cat, levels = c("High", "Moderate", "Low"))

# Arrange the data by API_cat and Tajima'sD (descending order for consistency)
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(API_cat) %>%
  arrange(desc(Pi), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Perform the linear regression
model <- lm(Pi ~ API, data = PF3D7_0711200_rowwise)
summary(model)

# Get the model coefficients
coefficients <- summary(model)$coefficients

slope <- coefficients["API", "Estimate"]
intercept <- coefficients["(Intercept)", "Estimate"]

# Create the formula text
formula_text <- paste("Pi = ", round(slope, 6), " * API + ", round(intercept, 6))

print(formula_text)

# Create the barplot for Pi with a line of best fit
ggplot(PF3D7_0711200_rowwise, aes(x = API, y = Pi, color = API_cat)) + 
  geom_point() +  # Add points to the plot
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a line of best fit
  labs(x = "API", y = "Pi", color = "API Category") +
  ggtitle("Global Nucleotide Diversity of PF3D7_0711200") +  # Add title here
  annotate("text", x = max(PF3D7_0711200_rowwise$API), y = max(PF3D7_0711200_rowwise$Pi), label = formula_text, 
           hjust = 1, vjust = 10, color = "black", size = 5) +  # Add formula text  
  theme_minimal() +
  scale_color_manual(values = c("High" = "darkred", "Moderate" = "darkorange", "Low" = "darkgreen")) +
  theme(
    plot.title = element_text(hjust = 0.5)  # Center the title
  )

##### Logarithmic Dot plot of Pi x API with line of best fit #####
library(broom)  # For tidy regression output

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat <- factor(PF3D7_0711200_rowwise$API_cat, levels = c("High", "Moderate", "Low"))

# Arrange the data by API_cat and Tajima'sD (descending order for consistency)
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(API_cat) %>%
  arrange(desc(Pi), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Perform the linear regression
model <- lm(Pi_pct ~ API, data = PF3D7_0711200_rowwise)
summary(model)

# Get the model coefficients
coefficients <- summary(model)$coefficients

slope <- coefficients["API", "Estimate"]
intercept <- coefficients["(Intercept)", "Estimate"]

# Create the formula text
formula_text <- paste("Pi (%) = ", round(slope, 6), " * API + ", round(intercept, 6))

print(formula_text)

# Create the barplot for Pi with a line of best fit and y axis normalized 
ggplot(PF3D7_0711200_rowwise, aes(x = API, y = Pi_pct, shape = Region)) + 
  geom_point(aes(fill = API_cat), color = "black", size = 4, stroke = 1.5) +  # Shape 21 supports both fill and border color
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black") +  # Add a line of best fit
  scale_x_log10() +  # Apply a logarithmic scale to the x-axis
  labs(x = "Log(API)", y = "Pi (%)", fill = "API Category", shape = "Region") +
  ggtitle("Regression of Global Nucleotide Diversity of \nConserved Pf Protein, Unknown Function (PF3D7_0711200)") +  # Add title here
  annotate("text", x = max(PF3D7_0711200_rowwise$API * 0.9), y = max(PF3D7_0711200_rowwise$Pi_pct * 2), label = formula_text, 
           hjust = 1, vjust = 0, color = "black", size = 5) +  # Add formula text  
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 0.65, by = .05), limits= c(0, 0.65)) +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +  # Fill color for API_cat
  scale_shape_manual(values = c("Central Africa" = 21,  # Circle
                                "East Africa" = 22,     # Triangle
                                "West Africa" = 23,     # Diamond
                                "South East Asia" = 24  # Square
  )) +  # Custom shapes for regions
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size   
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size =12) # Optional: tilt x-axis labels for readability
  )


# Create the barplot for Pi with a line of best fit
ggplot(PF3D7_0711200_rowwise, aes(x = API, y = Pi_pct, shape = Region)) + 
  geom_point(aes(fill = API_cat), color = "black", size = 4, stroke = 1.5) +  # Shape 21 supports both fill and border color
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black") +  # Add a line of best fit
  scale_x_log10() +  # Apply a logarithmic scale to the x-axis
  labs(x = "Log(API)", y = "Pi (%)", fill = "API Category", shape = "Region") +
  ggtitle("Regression of Global Nucleotide Diversity of Conserved Pf Protein, Unknown Function (PF3D7_0711200)") +  # Add title here
  annotate("text", x = max(PF3D7_0711200_rowwise$API), y = max(PF3D7_0711200_rowwise$Pi_pct), label = formula_text, 
           hjust = 1, vjust = 20, color = "black", size = 5) +  # Add formula text  
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +  # Fill color for API_cat
  scale_shape_manual(values = c("Central Africa" = 21,  # Circle
                                "East Africa" = 22,     # Triangle
                                "West Africa" = 23,     # Diamond
                                "South East Asia" = 24  # Square
  )) +  # Custom shapes for regions
  theme(
    plot.title = element_text(hjust = 0.5)  # Center the title
  )

# Create the scatterplot for Pi with a line of best fit and y axis normalized with colored API key
ggplot(PF3D7_0711200_rowwise, aes(x = API, y = Pi_pct, shape = Region, fill = API_cat)) + 
  geom_point(color = "black", size = 4, stroke = 1.5, show.legend = TRUE) +  # Shape 21 supports both fill and border color
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black") +  # Add a line of best fit
  scale_x_log10() +  # Apply a logarithmic scale to the x-axis
  labs(x = "Log(API)", y = "Pi (%)", fill = "API Category", shape = "Region") +
  ggtitle("Regression of Global Nucleotide Diversity of \nConserved Pf Protein, Unknown Function (PF3D7_0711200)") +  # Add title here
  annotate("text", x = max(PF3D7_0711200_rowwise$API * 0.9), y = max(PF3D7_0711200_rowwise$Pi_pct * 2), label = formula_text, 
           hjust = 1, vjust = 0, color = "black", size = 5) +  # Add formula text  
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 0.65, by = .05), limits= c(0, 0.65)) +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +  # Fill color for API_cat
  scale_shape_manual(values = c("Central Africa" = 21,  # Circle
                                "East Africa" = 22,     # Triangle
                                "West Africa" = 23,     # Diamond
                                "South East Asia" = 24  # Square
  )) +  # Custom shapes for regions
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size   
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size =12) # Optional: tilt x-axis labels for readability
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21, stroke = 0)))  # Remove stroke (border) in fill legend icons


##### Assessing Normality of Nucleotide/Haplotype Diversity ######

qqnorm(PF3D7_0711200_rowwise$Pi, main="Q-Q Plot of Average Nucleotide Diversity")
qqline(PF3D7_0711200_rowwise$Pi, col="darkred")

shapiro.test(PF3D7_0711200_rowwise$Pi)
shapiro.test(PF3D7_0711200_filtered$Pi)

qqnorm(PF3D7_0711200_rowwise$HaplotypeDiversity, main="Q-Q Plot of Average Haplotype Diversity")
qqline(PF3D7_0711200_rowwise$HaplotypeDiversity, col="darkred")

shapiro.test(PF3D7_0711200_rowwise$HaplotypeDiversity)
shapiro.test(PF3D7_0711200_filtered$HaplotypeDiversity)

##### Average Testing Nucleotide Diversity  ####

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/R_outputs/")
getwd()

# Perform ANOVA to test the difference in Pi across API_cat groups
anova_result_nuc <- aov(Pi ~ API_cat, data = PF3D7_0711200_rowwise)

# Summarize the ANOVA results
summary(anova_result_nuc)

anova_result_nuc_summary <- summary(anova_result_nuc)
anova_result_nuc_summary_table <- as.data.frame(anova_result_nuc_summary[[1]])
print(anova_result_nuc_summary_table)

#write.csv(anova_result_nuc_summary_table, "C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/R_outputs/PF3D7_0711200_anova_result_nuc_summary_table.csv", row.names = TRUE)

# Perform Tukey's HSD post-hoc test
tukey_result_nuc <- TukeyHSD(anova_result_nuc)
tukey_result_nuc_table <- as.data.frame(tukey_result_nuc$API_cat)  # Extract Tukey's HSD table
print(tukey_result_nuc_table)

#write.csv(tukey_result_nuc_table, "C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/R_outputs/PF3D7_0711200_tukey_result_nuc_table.csv", row.names = TRUE)

# View the post-hoc test results
tukey_result_nuc

##### Average Testing Nucleotide Diversity without Indonesia #####

# Filter out Indonesia from the dataset
PF3D7_0711200_filtered <- PF3D7_0711200_rowwise %>% filter(`Geographic.Site` != "Indonesia")

# Perform ANOVA to test the difference in Pi across API_cat groups, excluding Indonesia
anova_result_nuc_filtered <- aov(Pi ~ API_cat, data = PF3D7_0711200_filtered)

# Summarize the ANOVA results
summary(anova_result_nuc_filtered)

anova_result_nuc_filtered_summary <- summary(anova_result_nuc_filtered)
anova_result_nuc_filtered_summary_table <- as.data.frame(anova_result_nuc_filtered_summary[[1]])
print(anova_result_nuc_filtered_summary_table)

#write.csv(anova_result_nuc_filtered_summary_table, "C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/R_outputs/PF3D7_0711200_anova_result_nuc_filtered_summary_table.csv", row.names = TRUE)

# Perform Tukey's HSD post-hoc test
tukey_result_nuc_filtered <- TukeyHSD(anova_result_nuc_filtered)
tukey_result_nuc_filtered_table <- as.data.frame(tukey_result_nuc_filtered$API_cat)  # Extract Tukey's HSD table
print(tukey_result_nuc_filtered_table)

#write.csv(tukey_result_nuc_filtered_table, "C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/R_outputs/PF3D7_0711200_tukey_result_nuc_filtered_table.csv", row.names = TRUE)

# View the post-hoc test results
tukey_result_nuc_filtered

##### Average Testing Haplotype Diversity  ####

# Perform ANOVA to test the difference in HaplotypeDiversity across API_cat groups
anova_result_haplotype <- aov(HaplotypeDiversity ~ API_cat, data = PF3D7_0711200_rowwise)

# Summarize the ANOVA results
summary(anova_result_haplotype)

anova_result_haplotype_summary <- summary(anova_result_haplotype)
anova_result_haplotype_summary_table <- as.data.frame(anova_result_haplotype_summary[[1]])
print(anova_result_haplotype_summary_table)

#write.csv(anova_result_haplotype_summary_table, "C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/R_outputs/PF3D7_0711200_anova_result_haplotype_summary_table.csv", row.names = TRUE)

# Perform Tukey's HSD post-hoc test
tukey_result_haplotype <- TukeyHSD(anova_result_haplotype)
tukey_result_haplotype_table <- as.data.frame(tukey_result_haplotype$API_cat)  # Extract Tukey's HSD table
print(tukey_result_haplotype_table)

#write.csv(tukey_result_haplotype_table, "C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/R_outputs/PF3D7_0711200_tukey_result_haplotype_table.csv", row.names = TRUE)

# View the post-hoc test results
tukey_result_haplotype


##### Average Testing Haplotype Diversity without Indonesia #####

# Perform ANOVA to test the difference in HaplotypeDiversity across API_cat groups
anova_result_haplotype_filtered <- aov(HaplotypeDiversity ~ API_cat, data = PF3D7_0711200_filtered)

# Summarize the ANOVA results
summary(anova_result_haplotype_filtered)

anova_result_haplotype_filtered_summary <- summary(anova_result_haplotype_filtered)
anova_result_haplotype_filtered_summary_table <- as.data.frame(anova_result_haplotype_filtered_summary[[1]])
print(anova_result_haplotype_filtered_summary_table)

#write.csv(anova_result_haplotype_filtered_summary_table, "C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/R_outputs/PF3D7_0711200_anova_result_haplotype_filtered_summary_table.csv", row.names = TRUE)

# Perform Tukey's HSD post-hoc test
tukey_result_haplotype_filtered <- TukeyHSD(anova_result_haplotype_filtered)

tukey_result_haplotype_filtered_table <- as.data.frame(tukey_result_haplotype_filtered$API_cat)  # Extract Tukey's HSD table
print(tukey_result_haplotype_filtered_table)

#write.csv(tukey_result_haplotype_filtered_table, "C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/R_outputs/PF3D7_0711200_tukey_result_haplotype_filtered_table.csv", row.names = TRUE)

# View the post-hoc test results
tukey_result_haplotype_filtered

##### Median Testing Nucleotide Diversity ####
library(FSA)

kruskal.test(PF3D7_0711200_rowwise$Pi ~ PF3D7_0711200_rowwise$API_cat)

#Post-hoc testing for pairwise differences between medians 
dunnTest(Pi ~ API_cat, data=PF3D7_0711200_rowwise, method = "bonferroni")

##### Median Testing Nucleotide Diversity without Indonesia ####

kruskal.test(PF3D7_0711200_filtered$Pi ~ PF3D7_0711200_filtered$API_cat)

#Post-hoc testing for pairwise differences between medians 
dunnTest(Pi ~ API_cat, data=PF3D7_0711200_filtered, method = "bonferroni")

##### Median Testing Haplotype Diversity ####

kruskal.test(PF3D7_0711200_rowwise$HaplotypeDiversity ~ PF3D7_0711200_rowwise$API_cat)

#Post-hoc testing for pairwise differences between medians 
dunnTest(HaplotypeDiversity ~ API_cat, data=PF3D7_0711200_rowwise, method = "bonferroni")

##### Median Testing Haplotype Diversity without Indonesia ####

kruskal.test(PF3D7_0711200_filtered$HaplotypeDiversity ~ PF3D7_0711200_filtered$API_cat)

#Post-hoc testing for pairwise differences between medians 
dunnTest(HaplotypeDiversity ~ API_cat, data=PF3D7_0711200_filtered, method = "bonferroni")

