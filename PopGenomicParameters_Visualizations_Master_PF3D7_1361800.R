library(ggpattern) 
library(ggplot2)
library(dplyr)

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/R_inputs/")
getwd()

PF3D7_0711200_rowwise <- read.csv("PF3D7_1361800_popgenomic_estimates_rowwise.csv")

#### Removing country wide samples for countries with specific sites #####
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>% ###remove this if generating graphs for all sites
  filter(!Geographic.Site %in% c("Cambodia", "Vietnam", "Myanmar", "Senegal"))

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
                                       levels = c("West Africa", "Central Africa", "East Africa", "South East Asia"))

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
  ggtitle("Global Nucleotide Diversity of Glideosome-associated connector (PF3D7_1361800)") +  # Add title here
  scale_y_continuous(breaks = seq(0, 0.65, by = .05), limits= c(0, 0.65)) +
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +
  scale_pattern_manual(values = c("Central Africa" = "wave", 
                                  "East Africa" = "crosshatch", 
                                  "West Africa" = "stripe",  # You can set horizontal stripes here
                                  "South East Asia" = "circle")) + # or any other pattern
  theme(
    legend.position = "none",
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size 
    plot.title = element_text(hjust = 0.5),  # Center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Optional: tilt x-axis labels for readability
    axis.text.y = element_text(size =12)
  )

#### Nucleotide Diversity with Shading for Region and Adjusted API ############################

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat_adjusted <- factor(PF3D7_0711200_rowwise$API_cat_adjusted, levels = c("High", "Moderate", "Low"))

# Convert Region to a factor with the desired order
PF3D7_0711200_rowwise$Region <- factor(PF3D7_0711200_rowwise$Region, 
                                       levels = c("West Africa", "Central Africa", "East Africa", "South East Asia"))

# Arrange the data by Region, API_cat, and descending Pi
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(Region, API_cat_adjusted) %>%
  arrange(desc(Pi), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot with error bars and distinguish regions with patterns
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = Pi_pct, fill = API_cat_adjusted, pattern = Region)) + 
  geom_bar_pattern(stat = "identity", position = position_dodge(), width = 0.7, 
                   pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "white") + 
  geom_errorbar(aes(ymin = Pi_pct - PiSD_pct, ymax = Pi_pct + PiSD_pct), 
                width = 0.2, 
                position = position_dodge(0.7)) + 
  labs(x = "Geographic Site", y = "Pi (%)", fill = "API Category", pattern = "Region") +
  ggtitle("Global Nucleotide Diversity of Glideosome-associated connector (PF3D7_1361800)") +  # Add title here
  scale_y_continuous(breaks = seq(0, 0.65, by = .05), limits= c(0, 0.65)) +
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +
  scale_pattern_manual(values = c("Central Africa" = "wave", 
                                  "East Africa" = "crosshatch", 
                                  "West Africa" = "stripe",  # You can set horizontal stripes here
                                  "South East Asia" = "circle")) + # or any other pattern
  theme(
    legend.position = "none",
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size 
    plot.title = element_text(hjust = 0.5),  # Center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Optional: tilt x-axis labels for readability
    axis.text.y = element_text(size =12)
  )

#### Nucleotide Diversity with Shading for Region and Adjusted API and country in alphabetical order ############################

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat_adjusted <- factor(PF3D7_0711200_rowwise$API_cat_adjusted, levels = c("High", "Moderate", "Low"))

# Convert Region to a factor with the desired order
PF3D7_0711200_rowwise$Region <- factor(PF3D7_0711200_rowwise$Region, 
                                       levels = c("West Africa", "Central Africa", "East Africa", "South East Asia"))

#Change Geographic.site to character value so it can be sorted
PF3D7_0711200_rowwise$Geographic.Site <- as.character(PF3D7_0711200_rowwise$Geographic.Site)

# Arrange by API_cat_adjusted, Region, and Geographic.Site (alphabetically within groups)
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(API_cat_adjusted, Region) %>%
  arrange(API_cat_adjusted, Region, Geographic.Site) %>%
  ungroup()

# Convert Geographic.Site back to a factor after sorting
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, 
                                                levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot with error bars and distinguish regions with patterns
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = Pi_pct, fill = API_cat_adjusted, pattern = Region)) + 
  geom_bar_pattern(stat = "identity", position = position_dodge(), width = 0.7, 
                   pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "white") + 
  geom_errorbar(aes(ymin = Pi_pct - PiSD_pct, ymax = Pi_pct + PiSD_pct), 
                width = 0.2, 
                position = position_dodge(0.7)) + 
  labs(x = "Geographic Site", y = "Nucleotide diversity, π(%)", fill = "API Category", pattern = "Region") +
  #ggtitle("Global Nucleotide Diversity of Glideosome-associated connector (PF3D7_1361800)") +  # Add title here
  scale_y_continuous(breaks = seq(0, 0.65, by = .1), limits= c(0, 0.65)) +
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +
  scale_pattern_manual(values = c("Central Africa" = "wave", 
                                  "East Africa" = "crosshatch", 
                                  "West Africa" = "stripe",  # You can set horizontal stripes here
                                  "South East Asia" = "circle")) + # or any other pattern
  theme(
    legend.position = "none",
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size 
    plot.title = element_text(hjust = 0.5),  # Center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20), # Optional: tilt x-axis labels for readability
    axis.text.y = element_text(size =20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24, hjust = 1.25)
  )

##### Nucleotide Diversity Adjusted API and country in alphabetical order region shading removed #####

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat_adjusted <- factor(PF3D7_0711200_rowwise$API_cat_adjusted, levels = c("High", "Moderate", "Low"))

# Convert Region to a factor with the desired order
PF3D7_0711200_rowwise$Region <- factor(PF3D7_0711200_rowwise$Region, 
                                       levels = c("West Africa", "Central Africa", "East Africa", "South East Asia"))

# Change Geographic.site to character value so it can be sorted
PF3D7_0711200_rowwise$Geographic.Site <- as.character(PF3D7_0711200_rowwise$Geographic.Site)

# Arrange by API_cat_adjusted, Region, and Geographic.Site (alphabetically within groups)
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(API_cat_adjusted, Region) %>%
  arrange(API_cat_adjusted, Region, Geographic.Site) %>%
  ungroup()

# Convert Geographic.Site back to a factor after sorting
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, 
                                                levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot with error bars and color bars by API_cat_adjusted (no pattern)
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = Pi_pct, fill = API_cat_adjusted)) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +  # Remove the 'pattern' aesthetic
  geom_errorbar(aes(ymin = Pi_pct - PiSD_pct, ymax = Pi_pct + PiSD_pct), 
                width = 0.2, 
                position = position_dodge(0.7)) + 
  labs(x = "Geographic Site", y = "Nucleotide diversity, π(%)", fill = "API Category") +
  # ggtitle("Global Nucleotide Diversity of \nConserved Pf Protein, Unknown Function (PF3D7_0711200)") +  # Add title here
  scale_y_continuous(breaks = seq(0, 0.65, by = .1), limits = c(0, 0.65)) +
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +  # Color bars by API_cat_adjusted
  theme(
    legend.position = "none",
    legend.text = element_text(size = 13),  # Legend text size
    legend.title = element_text(size = 15),  # Legend title size 
    plot.title = element_text(hjust = 0.5),  # Center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),  # Optional: tilt x-axis labels for readability
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24, hjust = 1.25)
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
                                       levels = c("West Africa", "Central Africa", "East Africa", "South East Asia"))

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
  ggtitle("Global Haplotype Diversity of Glideosome-associated connector (PF3D7_1361800)") +  # Add title here
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +
  scale_pattern_manual(values = c("Central Africa" = "wave", 
                                  "East Africa" = "crosshatch", 
                                  "West Africa" = "stripe",  # You can set horizontal stripes here
                                  "South East Asia" = "circle")) + # or any other pattern
    scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1)) + # set y-axis limits and breaks
    theme(
    legend.position = "none",
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size   
    plot.title = element_text(hjust = 0.5), #center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size =12)) # Optional: tilt x-axis labels for readability

#### Haplotype Diversity with Shading for Region and Adjusted API ######################################

PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>% ###remove this if generating graphs for all sites
  filter(!Geographic.Site %in% c("Cambodia", "Vietnam", "Myanmar", "Senegal"))

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat_adjusted <- factor(PF3D7_0711200_rowwise$API_cat_adjusted, levels = c("High", "Moderate", "Low"))

# Convert Region to a factor with the desired order
PF3D7_0711200_rowwise$Region <- factor(PF3D7_0711200_rowwise$Region, 
                                       levels = c("West Africa", "Central Africa", "East Africa", "South East Asia"))

# Arrange the data by API_cat and descending Pi
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(Region, API_cat_adjusted) %>%
  arrange(desc(HaplotypeDiversity), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot with error bars
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = HaplotypeDiversity, fill = API_cat_adjusted, pattern = Region)) + 
  geom_bar_pattern(stat = "identity", position = position_dodge(), width = 0.7,
                   pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "white") + 
  geom_errorbar(aes(ymin = HaplotypeDiversity - HaplotypeSD, ymax = HaplotypeDiversity + HaplotypeSD), 
                width = 0.2, 
                position = position_dodge(0.7)) + 
  labs(x = "Geographic Site", y = "Haplotype Diversity", fill = "API Category", pattern = "Region") +
  ggtitle("Global Haplotype Diversity of Glideosome-associated connector (PF3D7_1361800)") +  # Add title here
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +
  scale_pattern_manual(values = c("Central Africa" = "wave", 
                                  "East Africa" = "crosshatch", 
                                  "West Africa" = "stripe",  # You can set horizontal stripes here
                                  "South East Asia" = "circle")) + # or any other pattern
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1)) + # set y-axis limits and breaks
  theme(
    legend.position = "none",
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size   
    plot.title = element_text(hjust = 0.5), #center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size =12)) # Optional: tilt x-axis labels for readability

#### Haplotype Diversity with Shading for Region and Adjusted API and country in alphabetical order ######################################

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat_adjusted <- factor(PF3D7_0711200_rowwise$API_cat_adjusted, levels = c("High", "Moderate", "Low"))

# Convert Region to a factor with the desired order
PF3D7_0711200_rowwise$Region <- factor(PF3D7_0711200_rowwise$Region, 
                                       levels = c("West Africa", "Central Africa", "East Africa", "South East Asia"))

#Change Geographic.site to character value so it can be sorted
PF3D7_0711200_rowwise$Geographic.Site <- as.character(PF3D7_0711200_rowwise$Geographic.Site)

# Arrange by API_cat_adjusted, Region, and Geographic.Site (alphabetically within groups)
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(API_cat_adjusted, Region) %>%
  arrange(API_cat_adjusted, Region, Geographic.Site) %>%
  ungroup()

# Convert Geographic.Site back to a factor after sorting
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, 
                                                levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot with error bars
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = HaplotypeDiversity, fill = API_cat_adjusted, pattern = Region)) + 
  geom_bar_pattern(stat = "identity", position = position_dodge(), width = 0.7,
                   pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "white") + 
  geom_errorbar(aes(ymin = HaplotypeDiversity - HaplotypeSD, ymax = HaplotypeDiversity + HaplotypeSD), 
                width = 0.2, 
                position = position_dodge(0.7)) + 
  labs(x = "Geographic Site", y = "Haplotype Diversity", fill = "API Category", pattern = "Region") +
  #ggtitle("Global Haplotype Diversity of Glideosome-associated connector (PF3D7_1361800)") +  # Add title here
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +
  scale_pattern_manual(values = c("Central Africa" = "wave", 
                                  "East Africa" = "crosshatch", 
                                  "West Africa" = "stripe",  # You can set horizontal stripes here
                                  "South East Asia" = "circle")) + # or any other pattern
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1)) + # set y-axis limits and breaks
  theme(
    legend.position = "none",
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size   
    plot.title = element_text(hjust = 0.5), #center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size =20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24)) # Optional: tilt x-axis labels for readability


##### Haplotype Diversity Adjusted API and country in alphabetical order region shading removed #####

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat_adjusted <- factor(PF3D7_0711200_rowwise$API_cat_adjusted, levels = c("High", "Moderate", "Low"))

# Convert Region to a factor with the desired order
PF3D7_0711200_rowwise$Region <- factor(PF3D7_0711200_rowwise$Region, 
                                       levels = c("West Africa", "Central Africa", "East Africa", "South East Asia"))

# Change Geographic.site to character value so it can be sorted
PF3D7_0711200_rowwise$Geographic.Site <- as.character(PF3D7_0711200_rowwise$Geographic.Site)

# Arrange by API_cat_adjusted, Region, and Geographic.Site (alphabetically within groups)
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(API_cat_adjusted, Region) %>%
  arrange(API_cat_adjusted, Region, Geographic.Site) %>%
  ungroup()

# Convert Geographic.Site back to a factor after sorting
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, 
                                                levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot with error bars
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = HaplotypeDiversity, fill = API_cat_adjusted)) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +  # Removed pattern aesthetic
  geom_errorbar(aes(ymin = HaplotypeDiversity - HaplotypeSD, ymax = HaplotypeDiversity + HaplotypeSD), 
                width = 0.2, 
                position = position_dodge(0.7)) + 
  labs(x = "Geographic Site", y = "Haplotype Diversity", fill = "API Category") + 
  # ggtitle("Global Haplotype Diversity of Glideosome-associated connector (PF3D7_1361800)") +  # Add title here
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +  # Color bars by API_cat_adjusted
  theme(
    legend.position = "none",
    legend.text = element_text(size = 13),  # Legend text size
    legend.title = element_text(size = 15),  # Legend title size   
    plot.title = element_text(hjust = 0.5),  # Center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24)  # Optional: tilt x-axis labels for readability
  )








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
                                       levels = c("West Africa", "Central Africa", "East Africa", "South East Asia"))

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
  geom_text(data = subset(PF3D7_0711200_rowwise, Tajima.sD > 2),
            aes(label = "*"), vjust = -0.5, color = "black", size = 8) +  # Place * above bars with Tajima.sD > 2
  geom_text(data = subset(PF3D7_0711200_rowwise, Tajima.sD < -2),
            aes(label = "*"), vjust = 1.5, color = "black", size = 8) +  # Place * below bars with Tajima.sD < -2
  labs(x = "Geographic Site", y = "Tajima's D", fill = "API Category", pattern = "Region") +
  ggtitle("Global Tajima's D of Glideosome-associated connector (PF3D7_1361800)") +  # Add title here
  scale_pattern_manual(values = c("Central Africa" = "wave", 
                                  "East Africa" = "crosshatch", 
                                  "West Africa" = "stripe",  # You can set horizontal stripes here
                                  "South East Asia" = "circle")) + # or any other pattern
  scale_y_continuous(breaks = seq(-3, 3, by = 0.5), limits= c(-3,3)) +
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size   
    plot.title = element_text(hjust = 0.5), #center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size =12) # Optional: tilt x-axis labels for readability
  )

#### Tajima's D with Region Shading and Adjusted API ############################################

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat_adjusted <- factor(PF3D7_0711200_rowwise$API_cat_adjusted, levels = c("High", "Moderate", "Low"))

# Convert Region to a factor with the desired order
PF3D7_0711200_rowwise$Region <- factor(PF3D7_0711200_rowwise$Region, 
                                       levels = c("West Africa", "Central Africa", "East Africa", "South East Asia"))

# Arrange the data by API_cat and Tajima'sD (descending order for consistency)
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(Region, API_cat_adjusted) %>%
  arrange(desc(Tajima.sD), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot for Tajima'sD
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = Tajima.sD, fill = API_cat_adjusted, pattern = Region)) + 
  geom_bar_pattern(stat = "identity", position = position_dodge(), width = 0.7,
                   pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "white") + 
  geom_hline(yintercept = 2, linetype = "dotted", color = "red", size = 1.4) +  # Add dotted red line at y = 2.5
  geom_hline(yintercept = -2, linetype = "dotted", color = "red", size = 1.4) + # Add dotted red line at y = -2.5
  geom_text(data = subset(PF3D7_0711200_rowwise, Tajima.sD > 2),
            aes(label = "+"), vjust = -0.5, color = "black", size = 8) +  # Place * above bars with Tajima.sD > 2
  geom_text(data = subset(PF3D7_0711200_rowwise, Tajima.sD < -2),
            aes(label = "+"), vjust = 1.5, color = "black", size = 8) +  # Place * below bars with Tajima.sD < -2
  labs(x = "Geographic Site", y = "Tajima's D", fill = "API Category", pattern = "Region") +
  ggtitle("Global Tajima's D of Glideosome-associated connector (PF3D7_1361800)") +  # Add title here
  scale_pattern_manual(values = c("Central Africa" = "wave", 
                                  "East Africa" = "crosshatch", 
                                  "West Africa" = "stripe",  # You can set horizontal stripes here
                                  "South East Asia" = "circle")) + # or any other pattern
  scale_y_continuous(breaks = seq(-3, 3, by = 0.5), limits= c(-3,3)) +
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size   
    plot.title = element_text(hjust = 0.5), #center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size =12) # Optional: tilt x-axis labels for readability
  )

#### Tajima's D with Region Shading and Adjusted API and country in alphabetical order ############################################

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat_adjusted <- factor(PF3D7_0711200_rowwise$API_cat_adjusted, levels = c("High", "Moderate", "Low"))

# Convert Region to a factor with the desired order
PF3D7_0711200_rowwise$Region <- factor(PF3D7_0711200_rowwise$Region, 
                                       levels = c("West Africa", "Central Africa", "East Africa", "South East Asia"))

#Change Geographic.site to character value so it can be sorted
PF3D7_0711200_rowwise$Geographic.Site <- as.character(PF3D7_0711200_rowwise$Geographic.Site)

# Arrange by API_cat_adjusted, Region, and Geographic.Site (alphabetically within groups)
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(API_cat_adjusted, Region) %>%
  arrange(API_cat_adjusted, Region, Geographic.Site) %>%
  ungroup()

# Convert Geographic.Site back to a factor after sorting
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, 
                                                levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot for Tajima'sD
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = Tajima.sD, fill = API_cat_adjusted, pattern = Region)) + 
  geom_bar_pattern(stat = "identity", position = position_dodge(), width = 0.7,
                   pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "white") + 
  geom_hline(yintercept = 2, linetype = "dotted", color = "red", size = 1.4) +  # Add dotted red line at y = 2.5
  geom_hline(yintercept = -2, linetype = "dotted", color = "red", size = 1.4) + # Add dotted red line at y = -2.5
  geom_text(data = subset(PF3D7_0711200_rowwise, Tajima.sD > 2),
            aes(label = "+"), vjust = -0.5, color = "black", size = 8) +  # Place * above bars with Tajima.sD > 2
  geom_text(data = subset(PF3D7_0711200_rowwise, Tajima.sD < -2),
            aes(label = "+"), vjust = 1.5, color = "black", size = 8) +  # Place * below bars with Tajima.sD < -2
  labs(x = "Geographic Site", y = "Tajima's D", fill = "API Category", pattern = "Region") +
  #ggtitle("Global Tajima's D of Glideosome-associated connector (PF3D7_1361800)") +  # Add title here
  scale_pattern_manual(values = c("Central Africa" = "wave", 
                                  "East Africa" = "crosshatch", 
                                  "West Africa" = "stripe",  # You can set horizontal stripes here
                                  "South East Asia" = "circle")) + # or any other pattern
  scale_y_continuous(breaks = seq(-3, 3, by = 1), limits= c(-3,3)) +
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size   
    plot.title = element_text(hjust = 0.5), #center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size =20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24) # Optional: tilt x-axis labels for readability
  )

##### Tajima's D  Adjusted API and country in alphabetical order region shading removed #####

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat_adjusted <- factor(PF3D7_0711200_rowwise$API_cat_adjusted, levels = c("High", "Moderate", "Low"))

# Convert Region to a factor with the desired order
PF3D7_0711200_rowwise$Region <- factor(PF3D7_0711200_rowwise$Region, 
                                       levels = c("West Africa", "Central Africa", "East Africa", "South East Asia"))

#Change Geographic.site to character value so it can be sorted
PF3D7_0711200_rowwise$Geographic.Site <- as.character(PF3D7_0711200_rowwise$Geographic.Site)

# Arrange by API_cat_adjusted, Region, and Geographic.Site (alphabetically within groups)
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(API_cat_adjusted, Region) %>%
  arrange(API_cat_adjusted, Region, Geographic.Site) %>%
  ungroup()

# Convert Geographic.Site back to a factor after sorting
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, 
                                                levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create the barplot for Tajima'sD
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = Tajima.sD, fill = API_cat_adjusted)) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) + 
  geom_hline(yintercept = 2, linetype = "dotted", color = "red", size = 1.4) +  # Add dotted red line at y = 2
  geom_hline(yintercept = -2, linetype = "dotted", color = "red", size = 1.4) + # Add dotted red line at y = -2
  geom_text(data = subset(PF3D7_0711200_rowwise, Tajima.sD > 2),
            aes(label = "+"), vjust = -0.5, color = "black", size = 8) +  # Place * above bars with Tajima.sD > 2
  geom_text(data = subset(PF3D7_0711200_rowwise, Tajima.sD < -2),
            aes(label = "+"), vjust = 1.5, color = "black", size = 8) +  # Place * below bars with Tajima.sD < -2
  labs(x = "Geographic Site", y = "Tajima's D", fill = "API Category") +
  scale_y_continuous(breaks = seq(-3, 3, by = 1), limits= c(-3,3)) +
  theme_minimal() +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 13), # Legend text size
    legend.title = element_text(size = 15), # Legend title size   
    plot.title = element_text(hjust = 0.5), # Center the title 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24), 
    axis.title.y = element_text(size = 24) # Optional: tilt x-axis labels for readability
  )



#### McDonald-Kreitman - Neutrality Index #####

#remove sites with Neutrality Index = NA

PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>% filter(!is.na(NeutralityIndex))

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat_adjusted <- factor(PF3D7_0711200_rowwise$API_cat_adjusted, levels = c("High", "Moderate", "Low"))

# Convert Region to a factor with the desired order
PF3D7_0711200_rowwise$Region <- factor(PF3D7_0711200_rowwise$Region, 
                                       levels = c("West Africa", "Central Africa", "East Africa", "South East Asia"))

#Change Geographic.site to character value so it can be sorted
PF3D7_0711200_rowwise$Geographic.Site <- as.character(PF3D7_0711200_rowwise$Geographic.Site)

# Arrange by API_cat_adjusted, Region, and Geographic.Site (alphabetically within groups)
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(API_cat_adjusted, Region) %>%
  arrange(API_cat_adjusted, Region, Geographic.Site) %>%
  ungroup()

# Convert Geographic.Site back to a factor after sorting
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, 
                                                levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Create plot for Neutrality Index with dots
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = NeutralityIndex, color = API_cat_adjusted, shape = Region)) + 
  geom_point(size = 4, position = position_jitter(width = 0.2, height = 0)) +  # Jitter to avoid overlap
  geom_hline(yintercept = 1, linetype = "dotted", color = "red", size = 1.4) +  
  labs(x = "Geographic Site", y = "Neutrality Index", color = "API Category", shape = "Region") +
  #ggtitle("Global Neutrality Index of of \nGlideosome-associated connector (PF3D7_1361800)") +  
  scale_y_continuous(breaks = seq(0, 15, by = 2.5), limits = c(0, 15)) +  
  theme_minimal() +
  scale_color_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray")) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),   
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size =20),
    axis.title.x = element_text(size =24),
    axis.title.y = element_text(size = 24)
  )

# MK with pseudo-log scale above 1 

# Load required package
library(scales)

# Function to create pseudo-log transformation
pseudo_log_trans <- function(base = 10) {
  trans_new("pseudo_log",
            transform = function(x) ifelse(x < 1, x, log10(x) + 1),  # Log scaling above 1
            inverse = function(x) ifelse(x < 1, x, base^(x - 1)),   # Reverse transformation
            breaks = c(0, 1, 10, 100),  # Break points at 0, 1, 10^1, 10^2
            format = function(x) ifelse(x < 1, as.character(x), parse(text = paste0("10^", log10(x))))
  )
}

# Create plot for Neutrality Index with pseudo-log y-axis
ggplot(PF3D7_0711200_rowwise, aes(x = Geographic.Site, y = NeutralityIndex, color = API_cat_adjusted, shape = Region)) + 
  geom_point(size = 4, position = position_jitter(width = 0.2, height = 0)) +  # Jitter to avoid overlap
  geom_hline(yintercept = 1, linetype = "dotted", color = "red", size = 1.4) +  
  labs(x = "Geographic Site", y = "Neutrality Index", color = "API Category", shape = "Region") +
  ggtitle("Global Neutrality Index of \nGlideosome-associated connector (PF3D7_1361800)") +  
  scale_y_continuous(
    trans = pseudo_log_trans(),  # Apply the pseudo-log transformation
    breaks = c(0, 1, 10, 100),   # Major breaks
    labels = c("0", "1", expression(10^1), expression(10^2)),  # Log-style labels
    limits = c(0, 100)  # Set y-axis range
  ) +  
  theme_minimal() +
  scale_color_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray")) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),   
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.minor = element_blank()  # Remove minor grid lines
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
#ggplot(PF3D7_0711200_rowwise, aes(x = API, y = Pi_pct, shape = Region)) + 
#  geom_point(aes(fill = API_cat), color = "black", size = 4, stroke = 1.5) +  # Shape 21 supports both fill and border color
#  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black") +  # Add a line of best fit
#  scale_x_log10() +  # Apply a logarithmic scale to the x-axis
#  labs(x = "Log(API)", y = "Pi (%)", fill = "API Category", shape = "Region") +
#  ggtitle("Regression of Global Nucleotide Diversity of \nConserved Pf Protein, Unknown Function (PF3D7_0711200)") +  # Add title here
#  annotate("text", x = max(PF3D7_0711200_rowwise$API * 0.9), y = max(PF3D7_0711200_rowwise$Pi_pct * 2), label = formula_text, 
#           hjust = 1, vjust = 0, color = "black", size = 5) +  # Add formula text  
#  theme_minimal() +
#  scale_y_continuous(breaks = seq(0, 0.65, by = .05), limits= c(0, 0.65)) +
#  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +  # Fill color for API_cat
#  scale_shape_manual(values = c("Central Africa" = 21,  # Circle
#                                "East Africa" = 22,     # Triangle
#                                "West Africa" = 23,     # Diamond
#                                "South East Asia" = 24  # Square
#  )) +  # Custom shapes for regions
#  theme(
#    plot.title = element_text(hjust = 0.5),  # Center the title
#    legend.text = element_text(size = 13), #legend text size
#    legend.title = element_text(size = 15), #legend title size   
#    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#    axis.text.y = element_text(size =12) # Optional: tilt x-axis labels for readability
#  )


# Create the barplot for Pi with a line of best fit
#ggplot(PF3D7_0711200_rowwise, aes(x = API, y = Pi_pct, shape = Region)) + 
#  geom_point(aes(fill = API_cat), color = "black", size = 4, stroke = 1.5) +  # Shape 21 supports both fill and border color
#  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black") +  # Add a line of best fit
#  scale_x_log10() +  # Apply a logarithmic scale to the x-axis
#  labs(x = "Log(API)", y = "Pi (%)", fill = "API Category", shape = "Region") +
#  ggtitle("Regression of Global Nucleotide Diversity of Conserved Pf Protein, Unknown Function (PF3D7_0711200)") +  # Add title here
#  annotate("text", x = max(PF3D7_0711200_rowwise$API), y = max(PF3D7_0711200_rowwise$Pi_pct), label = formula_text, 
#           hjust = 1, vjust = 20, color = "black", size = 5) +  # Add formula text  
#  theme_minimal() +
#  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +  # Fill color for API_cat
#  scale_shape_manual(values = c("Central Africa" = 21,  # Circle
#                                "East Africa" = 22,     # Triangle
#                                "West Africa" = 23,     # Diamond
#                                "South East Asia" = 24  # Square
#  )) +  # Custom shapes for regions
#  theme(
#    plot.title = element_text(hjust = 0.5)  # Center the title
#  )

# Create the scatterplot for Pi with a line of best fit and y axis normalized with colored API key
ggplot(PF3D7_0711200_rowwise, aes(x = API, y = Pi_pct, shape = Region, fill = API_cat)) + 
  geom_point(color = "black", size = 4, stroke = 1.5, show.legend = TRUE) +  # Shape 21 supports both fill and border color
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black") +  # Add a line of best fit
  scale_x_log10() +  # Apply a logarithmic scale to the x-axis
  labs(x = "Log(API)", y = "Pi (%)", fill = "API Category", shape = "Region") +
  ggtitle("Regression of Global Nucleotide Diversity of Glideosome-associated connector (PF3D7_1361800)") +  # Add title here
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
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),  # Center the title
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size   
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size =12) # Optional: tilt x-axis labels for readability
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21, stroke = 0)))  # Remove stroke (border) in fill legend icons

#### API Adjusted Regression #############################################

library(broom)  # For tidy regression output

# Convert API_cat to a factor with the desired order
PF3D7_0711200_rowwise$API_cat_adjusted <- factor(PF3D7_0711200_rowwise$API_cat_adjusted, levels = c("High", "Moderate", "Low"))

# Arrange the data by API_cat and Tajima'sD (descending order for consistency)
PF3D7_0711200_rowwise <- PF3D7_0711200_rowwise %>%
  group_by(API_cat_adjusted) %>%
  arrange(desc(Pi), .by_group = TRUE)

# Convert Geographic.Site to a factor with the correct order
PF3D7_0711200_rowwise$Geographic.Site <- factor(PF3D7_0711200_rowwise$Geographic.Site, levels = unique(PF3D7_0711200_rowwise$Geographic.Site))

# Perform the linear regression
model_API_adjusted <- lm(Pi_pct ~ API_regional, data = PF3D7_0711200_rowwise)
summary(model_API_adjusted)

# Get the model coefficients
coefficients_API_adjusted <- summary(model_API_adjusted)$coefficients

slope_API_adjusted <- coefficients_API_adjusted["API_regional", "Estimate"]
intercept_API_adjusted <- coefficients_API_adjusted["(Intercept)", "Estimate"]

# Create the formula text
formula_text_API_adjusted <- paste("π (%) = ", round(slope_API_adjusted, 6), " * API + ", round(intercept_API_adjusted, 6))

print(formula_text_API_adjusted)

# Create the scatterplot for Pi with a line of best fit and y axis normalized with colored API key
ggplot(PF3D7_0711200_rowwise, aes(x = API_regional, y = Pi_pct, shape = Region, fill = API_cat_adjusted)) + 
  geom_point(color = "black", size = 4, stroke = 1.5, show.legend = TRUE) +  # Shape 21 supports both fill and border color
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black") +  # Add a line of best fit
  scale_x_log10() +  # Apply a logarithmic scale to the x-axis
  labs(x = "Log(API)", y = "Nucleotide diversity, π(%)", fill = "API Category", shape = "Region") +
  #ggtitle("Regression of Global Nucleotide Diversity of \nGlideosome-associated connector (PF3D7_1361800)") +  # Add title here
  annotate("text", x = max(PF3D7_0711200_rowwise$API_regional * 0.9), y = max(PF3D7_0711200_rowwise$Pi_pct * 2), label = formula_text_API_adjusted, 
           hjust = 1, vjust = -1, color = "black", size = 7) +  # Add formula text  
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 0.65, by = .1), limits= c(0, 0.65)) +
  scale_fill_manual(values = c("High" = "deeppink3", "Moderate" = "darkorange", "Low" = "darkslategray1")) +  # Fill color for API_cat
  scale_shape_manual(values = c("Central Africa" = 21,  # Circle
                                "East Africa" = 22,     # Triangle
                                "West Africa" = 23,     # Diamond
                                "South East Asia" = 24  # Square
  )) +  # Custom shapes for regions
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),  # Center the title
    legend.text = element_text(size = 13), #legend text size
    legend.title = element_text(size = 15), #legend title size   
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size =20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24, )# Optional: tilt x-axis labels for readability
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21, stroke = 0)))  # Remove stroke (border) in fill legend icons







#### Calculcate R^2 value for original model #####
summary_model <- summary(model)
r_squared <- summary_model$r.squared
print(paste("R² =", round(r_squared, 4)))

# Extract the p-value for the API coefficient
p_value_non_log <- summary_model$coefficients["API", "Pr(>|t|)"]
print(paste("p-value for non-log-transformed API coefficient:", round(p_value_non_log, 6)))

#### Log transform API and calculcate R^2 ####

# Log-transform the API for the regression model
PF3D7_0711200_rowwise$log_API <- log10(PF3D7_0711200_rowwise$API)

# Perform the linear regression on the log-transformed API
log_model <- lm(Pi_pct ~ log_API, data = PF3D7_0711200_rowwise)
summary_log_model <- summary(log_model)

# Extract the R^2 and print it in the terminal
r_squared_log <- summary_log_model$r.squared
print(paste("R² =", round(r_squared_log, 4)))

#### Print p-value from model ####
p_value <- summary_log_model$coefficients["log_API", "Pr(>|t|)"]
print(paste("p-value for log_API coefficient:", round(p_value, 6)))

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




##### Adjusted API Median testing Nucleotide Diversity  #####

kruskal.test(PF3D7_0711200_rowwise$Pi ~ PF3D7_0711200_rowwise$API_cat_adjusted)

#Post-hoc testing for pairwise differences between medians 
dunnTest(Pi ~ API_cat_adjusted, data=PF3D7_0711200_rowwise, method = "bonferroni")

###### Adjusted API Median testing Haplotype Diversity ######

kruskal.test(PF3D7_0711200_rowwise$HaplotypeDiversity ~ PF3D7_0711200_rowwise$API_cat_adjusted)

#Post-hoc testing for pairwise differences between medians 
dunnTest(HaplotypeDiversity ~ API_cat_adjusted, data=PF3D7_0711200_rowwise, method = "bonferroni")



#Filter out Indonesia
PF3D7_0711200_filtered <- PF3D7_0711200_rowwise %>% filter(`Geographic.Site` != "Indonesia")

##### Median Testing Nucleotide Diversity without Indonesia ####
kruskal.test(PF3D7_0711200_filtered$Pi ~ PF3D7_0711200_filtered$API_cat_adjusted)

#Post-hoc testing for pairwise differences between medians 
dunnTest(Pi ~ API_cat_adjusted, data=PF3D7_0711200_filtered, method = "bonferroni")

###### Adjusted API Median testing Haplotype Diversity ######

kruskal.test(PF3D7_0711200_filtered$HaplotypeDiversity ~ PF3D7_0711200_filtered$API_cat_adjusted)

#Post-hoc testing for pairwise differences between medians 
dunnTest(HaplotypeDiversity ~ API_cat_adjusted, data=PF3D7_0711200_filtered, method = "bonferroni")

##### calculating group means and mean differences Nucleotide Diversity #####

# Calculate the mean Pi value for each API_cat_adjusted group
average_pi <- PF3D7_0711200_rowwise %>%
  group_by(API_cat_adjusted) %>%
  summarize(mean_Pi = mean(Pi, na.rm = TRUE))

# Print the results to a table
write.table(average_pi, "average_pi_by_group.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Create a table of mean differences between specified pairs
mean_differences_nuc <- tibble(
  Comparison = c("High - Moderate", "High - Low", "Moderate - Low"),
  Mean_Difference = c(
    average_pi$mean_Pi[average_pi$API_cat_adjusted == "High"] - average_pi$mean_Pi[average_pi$API_cat_adjusted == "Moderate"],
    average_pi$mean_Pi[average_pi$API_cat_adjusted == "High"] - average_pi$mean_Pi[average_pi$API_cat_adjusted == "Low"],
    average_pi$mean_Pi[average_pi$API_cat_adjusted == "Moderate"] - average_pi$mean_Pi[average_pi$API_cat_adjusted == "Low"]
  )
)

#Print the results to a table
write.table(mean_differences_nuc, "mean_differences_by_group.txt", sep = "\t", row.names = FALSE, quote = FALSE)


##### calculating group means and mean differences Haplotype Diversity #####

# Calculate the mean Pi value for each API_cat_adjusted group
average_hap <- PF3D7_0711200_rowwise %>%
  group_by(API_cat_adjusted) %>%
  summarize(mean_hap = mean(HaplotypeDiversity, na.rm = TRUE))

# Print the results to a table
write.table(average_hap, "average_haplotypediversity_by_group.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Create a table of mean differences between specified pairs
mean_differences_hap <- tibble(
  Comparison = c("High - Moderate", "High - Low", "Moderate - Low"),
  Mean_Difference = c(
    average_hap$mean_hap[average_hap$API_cat_adjusted == "High"] - average_hap$mean_hap[average_hap$API_cat_adjusted == "Moderate"],
    average_hap$mean_hap[average_hap$API_cat_adjusted == "High"] - average_hap$mean_hap[average_hap$API_cat_adjusted == "Low"],
    average_hap$mean_hap[average_hap$API_cat_adjusted == "Moderate"] - average_hap$mean_hap[average_hap$API_cat_adjusted == "Low"]
  )
)

# Print the results to a table
write.table(mean_differences_hap, "mean_differences_by_group.txt", sep = "\t", row.names = FALSE, quote = FALSE)