library(ggplot2)
library(dplyr)

###### PF3D7_0711200 #####

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/")
getwd()

Sliding_window_pi_df <- read.csv("PF3D7_0711200_sliding_window_pi_regions_window50_step25.csv")

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/nucleotide_index/")

epitope_count_df <- read.csv("PF3D7_0711200_22HLA_netmhcpan_predictions_nucleotide_index.csv") ### modify this as needed for inputs

# Remove any non-numeric characters (including leading/trailing spaces, tabs, etc.)
Sliding_window_pi_df$Midpoint <- gsub("[^0-9]", "", Sliding_window_pi_df$Midpoint)

# Convert the Midpoint column to numeric
Sliding_window_pi_df$Midpoint <- as.numeric(Sliding_window_pi_df$Midpoint)

#filter epitope count for positions greater than 0 
epitope_count_df <- epitope_count_df %>%
  filter(Epitope_count > 0) %>%
  mutate(y_min = 0, y_max = 0.03)  # Set y-coordinates for the colored bar


# Create the plot
ggplot(Sliding_window_pi_df, aes(x = Midpoint, y = Pi, group = Region, color = Region)) +
  geom_line() +            # Adds the line plot
  geom_segment(data = epitope_count_df, aes(x = Nucleotide_Position, xend = Nucleotide_Position, y = -0.001, yend = -0.0005),
               color = "orange", alpha = 0.5, size = 1, inherit.aes = FALSE) +
  theme_bw() +
  theme_minimal() +        # Optional: use a minimal theme
  labs(title = "Sliding Window Nucleotide Diversity in PF3D7_0711200",
       x = "Nucleotide Position", 
       y = "Nucleotide Diversity") + # Labels for x and y axes 
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels if necessary
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")) +
  scale_color_manual(values = c("West Africa" = "red2", "Central Africa" = "green2", "East Africa" = "black", 
                                "Southeast Asia" = "blue2"  )) +  # Optional: specify custom colors for each country
  scale_x_continuous(breaks = seq(0, max(Sliding_window_pi_df$Midpoint), by = 200)) + # Custom x-axis breaks
  scale_y_continuous(limits = c(-0.001, 0.03), breaks = seq(0.00, 0.03, by = 0.005)) #custom y -axis
  



#### PF3D7_0808100 #####

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/")
getwd()

Sliding_window_pi_df <- read.csv("PF3D7_0808100_sliding_window_pi_regions_window50_step25.csv")

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/nucleotide_index/")

epitope_count_df <- read.csv("PF3D7_0808100_22HLA_netmhcpan_predictions_nucleotide_index.csv") ### modify this as needed for inputs

# Remove any non-numeric characters (including leading/trailing spaces, tabs, etc.)
Sliding_window_pi_df$Midpoint <- gsub("[^0-9]", "", Sliding_window_pi_df$Midpoint)

# Convert the Midpoint column to numeric
Sliding_window_pi_df$Midpoint <- as.numeric(Sliding_window_pi_df$Midpoint)

#filter epitope count for positions greater than 0 
epitope_count_df <- epitope_count_df %>%
  filter(Epitope_count > 0) %>%
  mutate(y_min = 0, y_max = 0.03)  # Set y-coordinates for the colored bar


# Create the plot
ggplot(Sliding_window_pi_df, aes(x = Midpoint, y = Pi, group = Region, color = Region)) +
  geom_line() +            # Adds the line plot
  geom_segment(data = epitope_count_df, aes(x = Nucleotide_Position, xend = Nucleotide_Position, y = -0.001, yend = -0.0005),
               color = "orange", alpha = 0.5, size = 1, inherit.aes = FALSE) +
  theme_bw() +
  theme_minimal() +        # Optional: use a minimal theme
  labs(title = "Sliding Window Nucleotide Diversity in PF3D7_0808100",
       x = "Nucleotide Position", 
       y = "Nucleotide Diversity") + # Labels for x and y axes 
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels if necessary
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")) +
  scale_color_manual(values = c("West Africa" = "red2", "Central Africa" = "green2", "East Africa" = "black", 
                                "Southeast Asia" = "blue2"  )) +  # Optional: specify custom colors for each country
  scale_x_continuous(breaks = seq(0, max(Sliding_window_pi_df$Midpoint), by = 200)) + # Custom x-axis breaks
  scale_y_continuous(limits = c(-0.001, 0.03), breaks = seq(0.00, 0.03, by = 0.005)) #custom y -axis





#### PF3D7_1324300 #####

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/")
getwd()

Sliding_window_pi_df <- read.csv("PF3D7_1324300_sliding_window_pi_regions_window50_step25.csv")

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/nucleotide_index/")

epitope_count_df <- read.csv("PF3D7_1324300_22HLA_netmhcpan_predictions_nucleotide_index.csv") ### modify this as needed for inputs

# Remove any non-numeric characters (including leading/trailing spaces, tabs, etc.)
Sliding_window_pi_df$Midpoint <- gsub("[^0-9]", "", Sliding_window_pi_df$Midpoint)

# Convert the Midpoint column to numeric
Sliding_window_pi_df$Midpoint <- as.numeric(Sliding_window_pi_df$Midpoint)

#filter epitope count for positions greater than 0 
epitope_count_df <- epitope_count_df %>%
  filter(Epitope_count > 0) %>%
  mutate(y_min = 0, y_max = 0.03)  # Set y-coordinates for the colored bar


# Create the plot
ggplot(Sliding_window_pi_df, aes(x = Midpoint, y = Pi, group = Region, color = Region)) +
  geom_line() +            # Adds the line plot
  geom_segment(data = epitope_count_df, aes(x = Nucleotide_Position, xend = Nucleotide_Position, y = -0.001, yend = -0.0005),
               color = "orange", alpha = 0.5, size = 1, inherit.aes = FALSE) +
  theme_bw() +
  theme_minimal() +        # Optional: use a minimal theme
  labs(title = "Sliding Window Nucleotide Diversity in PF3D7_1324300",
       x = "Nucleotide Position", 
       y = "Nucleotide Diversity") + # Labels for x and y axes 
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels if necessary
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")) +
  scale_color_manual(values = c("West Africa" = "red2", "Central Africa" = "green2", "East Africa" = "black", 
                                "Southeast Asia" = "blue2"  )) +  # Optional: specify custom colors for each country
  scale_x_continuous(breaks = seq(0, max(Sliding_window_pi_df$Midpoint), by = 200)) + # Custom x-axis breaks
  scale_y_continuous(limits = c(-0.001, 0.03), breaks = seq(0.00, 0.03, by = 0.005)) #custom y -axis




#### PF3D7_1335900 #####

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/")
getwd()

Sliding_window_pi_df <- read.csv("PF3D7_1335900_sliding_window_pi_regions_window50_step25.csv")

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/nucleotide_index/")

epitope_count_df <- read.csv("PF3D7_1335900_22HLA_netmhcpan_predictions_nucleotide_index.csv") ### modify this as needed for inputs

# Remove any non-numeric characters (including leading/trailing spaces, tabs, etc.)
Sliding_window_pi_df$Midpoint <- gsub("[^0-9]", "", Sliding_window_pi_df$Midpoint)

# Convert the Midpoint column to numeric
Sliding_window_pi_df$Midpoint <- as.numeric(Sliding_window_pi_df$Midpoint)

#filter epitope count for positions greater than 0 
epitope_count_df <- epitope_count_df %>%
  filter(Epitope_count > 0) %>%
  mutate(y_min = 0, y_max = 0.03)  # Set y-coordinates for the colored bar


# Create the plot
ggplot(Sliding_window_pi_df, aes(x = Midpoint, y = Pi, group = Region, color = Region)) +
  geom_line() +            # Adds the line plot
  geom_segment(data = epitope_count_df, aes(x = Nucleotide_Position, xend = Nucleotide_Position, y = -0.001, yend = -0.0005),
               color = "orange", alpha = 0.5, size = 1, inherit.aes = FALSE) +
  theme_bw() +
  theme_minimal() +        # Optional: use a minimal theme
  labs(title = "Sliding Window Nucleotide Diversity in PF3D7_1335900",
       x = "Nucleotide Position", 
       y = "Nucleotide Diversity") + # Labels for x and y axes 
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels if necessary
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")) +
  scale_color_manual(values = c("West Africa" = "red2", "Central Africa" = "green2", "East Africa" = "black", 
                                "Southeast Asia" = "blue2"  )) +  # Optional: specify custom colors for each country
  scale_x_continuous(breaks = seq(0, max(Sliding_window_pi_df$Midpoint), by = 200)) + # Custom x-axis breaks
  scale_y_continuous(limits = c(-0.001, 0.03), breaks = seq(0.00, 0.03, by = 0.005)) #custom y -axis




#### PF3D7_1361800 #####

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/")
getwd()

Sliding_window_pi_df <- read.csv("PF3D7_1361800_sliding_window_pi_regions_window50_step25.csv")

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/nucleotide_index/")

epitope_count_df <- read.csv("PF3D7_1361800_22HLA_netmhcpan_predictions_nucleotide_index.csv") ### modify this as needed for inputs

# Remove any non-numeric characters (including leading/trailing spaces, tabs, etc.)
Sliding_window_pi_df$Midpoint <- gsub("[^0-9]", "", Sliding_window_pi_df$Midpoint)

# Convert the Midpoint column to numeric
Sliding_window_pi_df$Midpoint <- as.numeric(Sliding_window_pi_df$Midpoint)

#filter epitope count for positions greater than 0 
epitope_count_df <- epitope_count_df %>%
  filter(Epitope_count > 0) %>%
  mutate(y_min = 0, y_max = 0.03)  # Set y-coordinates for the colored bar


# Create the plot
ggplot(Sliding_window_pi_df, aes(x = Midpoint, y = Pi, group = Region, color = Region)) +
  geom_line() +            # Adds the line plot
  geom_segment(data = epitope_count_df, aes(x = Nucleotide_Position, xend = Nucleotide_Position, y = -0.001, yend = -0.0005),
               color = "orange", alpha = 0.5, size = 1, inherit.aes = FALSE) +
  theme_bw() +
  theme_minimal() +        # Optional: use a minimal theme
  labs(title = "Sliding Window Nucleotide Diversity in PF3D7_1361800",
       x = "Nucleotide Position", 
       y = "Nucleotide Diversity") + # Labels for x and y axes 
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels if necessary
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")) +
  scale_color_manual(values = c("West Africa" = "red2", "Central Africa" = "green2", "East Africa" = "black", 
                                "Southeast Asia" = "blue2"  )) +  # Optional: specify custom colors for each country
  scale_x_continuous(breaks = seq(0, max(Sliding_window_pi_df$Midpoint), by = 200)) + # Custom x-axis breaks
  scale_y_continuous(limits = c(-0.001, 0.03), breaks = seq(0.00, 0.03, by = 0.005)) #custom y -axis













####### Section with Epitope_count categories ######

###### PF3D7_0711200 #####

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/")
getwd()

Sliding_window_pi_df <- read.csv("PF3D7_0711200_sliding_window_pi_regions_window50_step25.csv")

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/nucleotide_index/16HLA_pf3k/")

epitope_count_df <- read.csv("PF3D7_0711200_16HLAforpf3k_netmhcpan_predictions_nucleotide_index.csv") ### modify this as needed for inputs

# Remove any non-numeric characters (including leading/trailing spaces, tabs, etc.)
Sliding_window_pi_df$Midpoint <- gsub("[^0-9]", "", Sliding_window_pi_df$Midpoint)

# Convert the Midpoint column to numeric
Sliding_window_pi_df$Midpoint <- as.numeric(Sliding_window_pi_df$Midpoint)

#filter epitope count for positions greater than 0 
epitope_count_df <- epitope_count_df %>%
  filter(Epitope_count > 0) %>%
  mutate(y_min = 0, y_max = 0.03)  # Set y-coordinates for the colored bar

# Add a new column for the color category based on Epitope_count ranges
epitope_count_df <- epitope_count_df %>%
  mutate(Epitope_color = case_when(
    Epitope_count >= 1 & Epitope_count <= 5 ~ "purple",
    Epitope_count >= 6 & Epitope_count <= 10 ~ "magenta",
    Epitope_count >= 11 & Epitope_count <= 16 ~ "darkturquoise",
    TRUE ~ NA_character_  # Handle any unexpected values
  ))


# Create the plot
ggplot(Sliding_window_pi_df, aes(x = Midpoint, y = Pi, group = Region)) +
  geom_line(aes(color = Region)) +            # Adds the line plot with Region-based colors
  geom_segment(data = epitope_count_df, 
               aes(x = Nucleotide_Position, xend = Nucleotide_Position, y = -0.001, yend = -0.0005, color = Epitope_color),
               alpha = 0.5, size = 1, inherit.aes = FALSE) +  # Bar colored by Epitope_count range
  theme_bw() +
  theme_minimal() +        # Optional: use a minimal theme
  #labs(title = "Sliding Window Nucleotide Diversity in PF3D7_0711200",
  labs(
       x = "Nucleotide Position", 
       y = "Nucleotide Diversity") + # Labels for x and y axes 
  theme(
    axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    legend.position = "none",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")) +
  scale_color_manual(
    values = c("West Africa" = "red2", "Central Africa" = "green2", 
               "East Africa" = "black", "Southeast Asia" = "blue2", 
               "purple" = "purple", "magenta" = "magenta", "darkturquoise" = "darkturquoise"),
    aesthetics = c("color"),  # Apply to both layers
    guide = guide_legend(title = "Region / Epitope Count")) +  # Combined legend
  scale_x_continuous(breaks = seq(0, max(Sliding_window_pi_df$Midpoint), by = 200)) + # Custom x-axis breaks
  scale_y_continuous(limits = c(-0.001, 0.03), breaks = seq(0.00, 0.03, by = 0.005)) # Custom y-axis





##### PF3D7_0808100 #####

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/")
getwd()

Sliding_window_pi_df <- read.csv("PF3D7_0808100_sliding_window_pi_regions_window50_step25.csv")

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/nucleotide_index/16HLA_pf3k/")

epitope_count_df <- read.csv("PF3D7_0808100_16HLAforpf3k_netmhcpan_predictions_nucleotide_index.csv") ### modify this as needed for inputs

# Remove any non-numeric characters (including leading/trailing spaces, tabs, etc.)
Sliding_window_pi_df$Midpoint <- gsub("[^0-9]", "", Sliding_window_pi_df$Midpoint)

# Convert the Midpoint column to numeric
Sliding_window_pi_df$Midpoint <- as.numeric(Sliding_window_pi_df$Midpoint)

#filter epitope count for positions greater than 0 
epitope_count_df <- epitope_count_df %>%
  filter(Epitope_count > 0) %>%
  mutate(y_min = 0, y_max = 0.03)  # Set y-coordinates for the colored bar

# Add a new column for the color category based on Epitope_count ranges
epitope_count_df <- epitope_count_df %>%
  mutate(Epitope_color = case_when(
    Epitope_count >= 1 & Epitope_count <= 5 ~ "purple",
    Epitope_count >= 6 & Epitope_count <= 10 ~ "magenta",
    Epitope_count >= 11 & Epitope_count <= 16 ~ "darkturquoise",
    TRUE ~ NA_character_  # Handle any unexpected values
  ))


# Create the plot
ggplot(Sliding_window_pi_df, aes(x = Midpoint, y = Pi, group = Region)) +
  geom_line(aes(color = Region)) +            # Adds the line plot with Region-based colors
  geom_segment(data = epitope_count_df, 
               aes(x = Nucleotide_Position, xend = Nucleotide_Position, y = -0.001, yend = -0.0005, color = Epitope_color),
               alpha = 0.5, size = 1, inherit.aes = FALSE) +  # Bar colored by Epitope_count range
  theme_bw() +
  theme_minimal() +        # Optional: use a minimal theme
  #labs(title = "Sliding Window Nucleotide Diversity in PF3D7_0808100",
  labs(     
       x = "Nucleotide Position", 
       y = "Nucleotide Diversity") + # Labels for x and y axes 
  theme(
    axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    legend.position = "none",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")) +
  scale_color_manual(
    values = c("West Africa" = "red2", "Central Africa" = "green2", 
               "East Africa" = "black", "Southeast Asia" = "blue2", 
               "purple" = "purple", "magenta" = "magenta", "darkturquoise" = "darkturquoise"),
    aesthetics = c("color"),  # Apply to both layers
    guide = guide_legend(title = "Region / Epitope Count")) +  # Combined legend
  scale_x_continuous(breaks = seq(0, max(Sliding_window_pi_df$Midpoint), by = 200)) + # Custom x-axis breaks
  scale_y_continuous(limits = c(-0.001, 0.03), breaks = seq(0.00, 0.03, by = 0.005)) # Custom y-axis




##### PF3D7_1324300 #####

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/")
getwd()

Sliding_window_pi_df <- read.csv("PF3D7_1324300_sliding_window_pi_regions_window50_step25.csv")

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/nucleotide_index/16HLA_pf3k/")

epitope_count_df <- read.csv("PF3D7_1324300_16HLAforpf3k_netmhcpan_predictions_nucleotide_index.csv") ### modify this as needed for inputs

# Remove any non-numeric characters (including leading/trailing spaces, tabs, etc.)
Sliding_window_pi_df$Midpoint <- gsub("[^0-9]", "", Sliding_window_pi_df$Midpoint)

# Convert the Midpoint column to numeric
Sliding_window_pi_df$Midpoint <- as.numeric(Sliding_window_pi_df$Midpoint)

#filter epitope count for positions greater than 0 
epitope_count_df <- epitope_count_df %>%
  filter(Epitope_count > 0) %>%
  mutate(y_min = 0, y_max = 0.03)  # Set y-coordinates for the colored bar

# Add a new column for the color category based on Epitope_count ranges
epitope_count_df <- epitope_count_df %>%
  mutate(Epitope_color = case_when(
    Epitope_count >= 1 & Epitope_count <= 5 ~ "purple",
    Epitope_count >= 6 & Epitope_count <= 10 ~ "magenta",
    Epitope_count >= 11 & Epitope_count <= 16 ~ "darkturquoise",
    TRUE ~ NA_character_  # Handle any unexpected values
  ))


# Create the plot
ggplot(Sliding_window_pi_df, aes(x = Midpoint, y = Pi, group = Region)) +
  geom_line(aes(color = Region)) +            # Adds the line plot with Region-based colors
  geom_segment(data = epitope_count_df, 
               aes(x = Nucleotide_Position, xend = Nucleotide_Position, y = -0.001, yend = -0.0005, color = Epitope_color),
               alpha = 0.5, size = 1, inherit.aes = FALSE) +  # Bar colored by Epitope_count range
  theme_bw() +
  theme_minimal() +        # Optional: use a minimal theme
  #labs(title = "Sliding Window Nucleotide Diversity in PF3D7_1324300",
  labs( 
      x = "Nucleotide Position", 
       y = "Nucleotide Diversity") + # Labels for x and y axes 
  theme(
    axis.text.x = element_text(size = 20, angle = 90, hjust = 1),  
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    legend.position = "none",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")) +
  scale_color_manual(
    values = c("West Africa" = "red2", "Central Africa" = "green2", 
               "East Africa" = "black", "Southeast Asia" = "blue2", 
               "purple" = "purple", "magenta" = "magenta", "darkturquoise" = "darkturquoise"),
    aesthetics = c("color"),  # Apply to both layers
    guide = guide_legend(title = "Region / Epitope Count")) +  # Combined legend
  scale_x_continuous(breaks = seq(0, max(Sliding_window_pi_df$Midpoint), by = 1000)) + # Custom x-axis breaks
  scale_y_continuous(limits = c(-0.001, 0.03), breaks = seq(0.00, 0.03, by = 0.005)) # Custom y-axis



#### PF3D7_1335900 #####

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/")
getwd()

Sliding_window_pi_df <- read.csv("PF3D7_1335900_sliding_window_pi_regions_window50_step25.csv")

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/nucleotide_index/16HLA_pf3k/")

epitope_count_df <- read.csv("PF3D7_1335900_16HLAforpf3k_netmhcpan_predictions_nucleotide_index.csv") ### modify this as needed for inputs

# Remove any non-numeric characters (including leading/trailing spaces, tabs, etc.)
Sliding_window_pi_df$Midpoint <- gsub("[^0-9]", "", Sliding_window_pi_df$Midpoint)

# Convert the Midpoint column to numeric
Sliding_window_pi_df$Midpoint <- as.numeric(Sliding_window_pi_df$Midpoint)

#filter epitope count for positions greater than 0 
epitope_count_df <- epitope_count_df %>%
  filter(Epitope_count > 0) %>%
  mutate(y_min = 0, y_max = 0.03)  # Set y-coordinates for the colored bar

# Add a new column for the color category based on Epitope_count ranges
epitope_count_df <- epitope_count_df %>%
  mutate(Epitope_color = case_when(
    Epitope_count >= 1 & Epitope_count <= 5 ~ "purple",
    Epitope_count >= 6 & Epitope_count <= 10 ~ "magenta",
    Epitope_count >= 11 & Epitope_count <= 16 ~ "darkturquoise",
    TRUE ~ NA_character_  # Handle any unexpected values
  ))


# Create the plot
ggplot(Sliding_window_pi_df, aes(x = Midpoint, y = Pi, group = Region)) +
  geom_line(aes(color = Region)) +            # Adds the line plot with Region-based colors
  geom_segment(data = epitope_count_df, 
               aes(x = Nucleotide_Position, xend = Nucleotide_Position, y = -0.001, yend = -0.0005, color = Epitope_color),
               alpha = 0.5, size = 1, inherit.aes = FALSE) +  # Bar colored by Epitope_count range
  theme_bw() +
  theme_minimal() +        # Optional: use a minimal theme
  #labs(title = "Sliding Window Nucleotide Diversity in PF3D7_1335900",
  labs( 
      x = "Nucleotide Position", 
       y = "Nucleotide Diversity") + # Labels for x and y axes 
  theme(
    axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    legend.position = "none",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")) +
  scale_color_manual(
    values = c("West Africa" = "red2", "Central Africa" = "green2", 
               "East Africa" = "black", "Southeast Asia" = "blue2", 
               "purple" = "purple", "magenta" = "magenta", "darkturquoise" = "darkturquoise"),
    aesthetics = c("color"),  # Apply to both layers
    guide = guide_legend(title = "Region / Epitope Count")) +  # Combined legend
  scale_x_continuous(breaks = seq(0, max(Sliding_window_pi_df$Midpoint), by = 100)) + # Custom x-axis breaks
  scale_y_continuous(limits = c(-0.001, 0.03), breaks = seq(0.00, 0.03, by = 0.005)) # Custom y-axis






##### PF3D7_1361800 #####

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/")
getwd()

Sliding_window_pi_df <- read.csv("PF3D7_1361800_sliding_window_pi_regions_window50_step25.csv")

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/Sliding_window_Pi_analysis/nucleotide_index/16HLA_pf3k/")

epitope_count_df <- read.csv("PF3D7_1361800_16HLAforpf3k_netmhcpan_predictions_nucleotide_index.csv") ### modify this as needed for inputs

# Remove any non-numeric characters (including leading/trailing spaces, tabs, etc.)
Sliding_window_pi_df$Midpoint <- gsub("[^0-9]", "", Sliding_window_pi_df$Midpoint)

# Convert the Midpoint column to numeric
Sliding_window_pi_df$Midpoint <- as.numeric(Sliding_window_pi_df$Midpoint)

#filter epitope count for positions greater than 0 
epitope_count_df <- epitope_count_df %>%
  filter(Epitope_count > 0) %>%
  mutate(y_min = 0, y_max = 0.03)  # Set y-coordinates for the colored bar

# Add a new column for the color category based on Epitope_count ranges
epitope_count_df <- epitope_count_df %>%
  mutate(Epitope_color = case_when(
    Epitope_count >= 1 & Epitope_count <= 5 ~ "purple",
    Epitope_count >= 6 & Epitope_count <= 10 ~ "magenta",
    Epitope_count >= 11 & Epitope_count <= 16 ~ "darkturquoise",
    TRUE ~ NA_character_  # Handle any unexpected values
  ))


# Create the plot
ggplot(Sliding_window_pi_df, aes(x = Midpoint, y = Pi, group = Region)) +
  geom_line(aes(color = Region)) +            # Adds the line plot with Region-based colors
  geom_segment(data = epitope_count_df, 
               aes(x = Nucleotide_Position, xend = Nucleotide_Position, y = -0.001, yend = -0.0005, color = Epitope_color),
               alpha = 0.5, size = 1, inherit.aes = FALSE) +  # Bar colored by Epitope_count range
  theme_bw() +
  theme_minimal() +        # Optional: use a minimal theme
  #labs(title = "Sliding Window Nucleotide Diversity in PF3D7_1361800",
  labs( 
      x = "Nucleotide Position", 
       y = "Nucleotide Diversity") + # Labels for x and y axes 
  theme(
    axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size =24),
    axis.title.y = element_text(size = 24),
    legend.position = "none",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")) +
  scale_color_manual(
    values = c("West Africa" = "red2", "Central Africa" = "green2", 
               "East Africa" = "black", "Southeast Asia" = "blue2", 
               "purple" = "purple", "magenta" = "magenta", "darkturquoise" = "darkturquoise"),
    aesthetics = c("color"),  # Apply to both layers
    guide = guide_legend(title = "Region / Epitope Count")) +  # Combined legend
  scale_x_continuous(breaks = seq(0, max(Sliding_window_pi_df$Midpoint), by = 500)) + # Custom x-axis breaks
  scale_y_continuous(limits = c(-0.001, 0.03), breaks = seq(0.00, 0.03, by = 0.005)) # Custom y-axis



