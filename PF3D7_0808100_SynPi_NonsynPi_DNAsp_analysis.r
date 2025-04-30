library(ggplot2)
library(tidyr)

setwd("C:/Users/Ryan/Desktop/CarneiroDaSilva/Pf3k/PopGenomic_Estimates/DNAsp_based_assessment_of_Pi_in_NONSYN_SYN/")
getwd()

Locus_nonsyn_syn_df <- read.csv("PF3D7_0808100_Pi_nonsyn_syn_estimates_by_country_rowwise_reduced.csv")

# Create Pi_pct variables #
Locus_nonsyn_syn_df <- Locus_nonsyn_syn_df %>%
  mutate(
    Synonymous = Syn_pi * 100,
    Nonsynonymous = Nonsyn_pi * 100
  )

# Reshape the data to long format
Locus_nonsyn_syn_df_long <- Locus_nonsyn_syn_df %>%
  pivot_longer(
    cols = c(Synonymous, Nonsynonymous),
    names_to = "Type",
    values_to = "Value"
  )

# Convert type to factor variable 
Locus_nonsyn_syn_df_long$Type <- as.factor(Locus_nonsyn_syn_df_long$Type)

# Assign region colors
region_colors <- c("Southeast Asia" = "blue3", "West Africa" = "red3", "Central Africa" = "azure4", "East Africa" = "green3")

# Define the custom order of regions
region_order <- c("West Africa", "Central Africa", "East Africa", "Southeast Asia")

# Ensure the Region column follows the correct order
Locus_nonsyn_syn_df_long$Region <- factor(Locus_nonsyn_syn_df_long$Region, levels = region_order)

# Sort countries alphabetically within each region
Locus_nonsyn_syn_df_long$Country <- factor(
  Locus_nonsyn_syn_df_long$Country,
  levels = unique(Locus_nonsyn_syn_df_long$Country[order(Locus_nonsyn_syn_df_long$Region, Locus_nonsyn_syn_df_long$Country)]),
  ordered = TRUE
)



##### Create a double bar plot showing nonsyn and syn values for each country ######

# Create the bar plot with tilted and properly aligned country labels
ggplot(Locus_nonsyn_syn_df_long, aes(x = Country, y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(
    #title = "Nonsynonymous and Synonymous Nucleotide Diversity in PF3D7_0808100",
    x = "Country",
    y = "Nucleotide Diversity (%)",
    fill = "Position Type"
  ) +
  scale_fill_manual(values = c("Synonymous" = "aquamarine3", "Nonsynonymous" = "brown2")) +
  theme_minimal() +
  theme(
    legend.position = "none", # to bring back legend just delete this
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24, hjust = 1.2)
  ) +
  # Add region lines
  geom_segment(data = unique(Locus_nonsyn_syn_df_long[c("Country", "Region")]),
               aes(x = as.numeric(Country) - 0.5, 
                   xend = as.numeric(Country) + 0.5, 
                   y = -0.01, yend = -0.01, color = Region), ### 0.00001 puts bar touching bottom of bar plot bars
               size = 3, inherit.aes = FALSE) +
  scale_color_manual(values = region_colors, guide = "none") +
  scale_color_manual(values = region_colors, guide = "legend") +  # Add region color to the legend
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) + # Order the legends
  scale_y_continuous( #customize y-axis here
    limits = c(-.01, 0.6), #min/max
    breaks = seq(0, 1.0, by = 0.1) #set breaks
  )

### Create Nonsyn_pi / syn_pi variable

Locus_nonsyn_syn_df <- Locus_nonsyn_syn_df %>%
  mutate(Ratio_nonsyn_pi_syn_pi = Nonsyn_pi / Syn_pi)

mean(Locus_nonsyn_syn_df$Ratio_nonsyn_pi_syn_pi)