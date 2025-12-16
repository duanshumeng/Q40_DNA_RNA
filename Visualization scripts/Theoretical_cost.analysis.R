# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggpubr)
library(forcats)

setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/文章/GB投稿/Reviewer_comments_reply/Figures')
# Create data frame
seq_data <- data.frame(
  Application_scenarios = c("RNA-seq", "RNA-seq", 
                            "Germline Variant\nWES", "Germline Variant\nWES",
                            "Somatic Variant\nWES", "Somatic Variant\nWES",
                            "Germline Variant\nWGS", "Germline Variant\nWGS",
                            "Somatic Variant\nWGS", "Somatic Variant\nWGS"),
  Type = rep(c("Q30", "Q40"), 5),
  Volume_GB = c(10, 4, 1.05, 0.7, 3.15, 2.1, 96, 64, 288, 192),
  Total_Cost_RMB = c(600, 420, 331.5, 321, 394.5, 363, 3180, 2220, 8940, 6060),
  Data_volume_compression_ratio = c(NA, 60, NA, 33.3, NA, 33.3, NA, 33.3, NA, 33.3),
  Cost_compression_ratio = c(NA, 24, NA, 2.2, NA, 5.8, NA, 28.8, NA, 31.7)
)

# Create grouping variables
seq_data$Group <- factor(seq_data$Application_scenarios, 
                         levels = unique(seq_data$Application_scenarios))
seq_data$Quality <- factor(seq_data$Type, levels = c("Q30", "Q40"))

# Compute change percentages
seq_data <- seq_data %>%
  dplyr::group_by(Application_scenarios) %>%
  dplyr::mutate(
    Volume_change_pct = (Volume_GB[Quality == "Q30"] - Volume_GB[Quality == "Q40"]) / Volume_GB[Quality == "Q30"] * 100,
    Cost_change_pct = (Total_Cost_RMB[Quality == "Q30"] - Total_Cost_RMB[Quality == "Q40"]) / Total_Cost_RMB[Quality == "Q30"] * 100
  ) %>%
  ungroup()

# Nature-style theme settings
nature_theme <- theme_minimal(base_size = 11) +
  theme(
    text = element_text(color = "#333333"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, 
                              margin = margin(b = 15)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#666666",
                                 margin = margin(b = 20)),
    axis.title = element_text(size = 11, face = "plain"),
    axis.text = element_text(size = 10, color = "#333333"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.line = element_line(color = "#333333", size = 0.5),
    panel.grid.major = element_line(color = "#E0E0E0", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.position = "top",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(20, 20, 20, 20),
    strip.text = element_text(size = 10, face = "bold")
  )


# 3. Visualize change percentages (heatmap style)
change_data <- seq_data %>%
  filter(Quality == "Q40") %>%
  mutate(
    Data_volume_change = -Volume_change_pct,  # convert to positive to indicate compression
    Cost_saving = -Cost_change_pct            # convert to positive to indicate savings
  ) %>%
  select(Application_scenarios, Data_volume_change, Cost_saving) %>%
  pivot_longer(cols = c(Data_volume_change, Cost_saving), 
               names_to = "Metric", values_to = "Value")

p3 <- ggplot(change_data, aes(x = Application_scenarios, y = Metric, fill = Value)) +
  geom_tile(color = "white", size = 1.5) +
  geom_text(aes(label = paste0(round(Value, 1), "%")), 
            color = "black", size = 4.5, fontface = "bold") +
  scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", 
                       midpoint = -10, 
                       limits = c(-60, 0),
                       name = "Change (%)",
                       breaks = c(-60,-50,-40,-30,-20,-10,0)) +
  scale_y_discrete(labels = c("Data Volume\nCompression", "Cost\nSavings")) +
  labs(title = "Impact of Quality Improvement: Q30 to Q40",
       subtitle = "",
       x = "Application Scenario",
       y = "Metric") +
  nature_theme +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9))

ggsave('Extend-Fig6.cost-saving.pdf',p3,width=8.15, height=5.7,dpi = 300)