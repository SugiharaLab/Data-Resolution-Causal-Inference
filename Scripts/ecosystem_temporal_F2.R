library(ggplot2)
library(gridExtra)

# Define a function to load data and identify resolved interactions
identify_resolved_interactions <- function(filepath) {
  load(filepath)
  resolved <- ccm_output_combined$convergence_p_val < 0.05 & ccm_output_combined$ccm_rho > ccm_output_combined$`0.95P_trials`
  ccm_output_combined[resolved, c("target", "lib")]
}

# Systems to analyze
systems <- c("PEB", "LZ", "NS")
results <- data.frame(system = character(), type = character(), count = integer())

# Loop through each system
for (s in systems) {
  monthly_path <- paste0(getwd(), "/Data/CCM Output/", s, "_CCM_monthly.RData")
  annual_path <- paste0(getwd(), "/Data/CCM Output/", s, "_CCM_annual.RData")
  
  monthly_resolved <- identify_resolved_interactions(monthly_path)
  annual_resolved <- identify_resolved_interactions(annual_path)
  
  exclusive_monthly <- setdiff(monthly_resolved, annual_resolved)
  exclusive_annual <- setdiff(annual_resolved, monthly_resolved)
  common_resolved <- intersect(monthly_resolved, annual_resolved)
  
  results <- rbind(results, data.frame(system = s, type = "Monthly Only", count = nrow(exclusive_monthly)))
  results <- rbind(results, data.frame(system = s, type = "Annual Only", count = nrow(exclusive_annual)))
  results <- rbind(results, data.frame(system = s, type = "Both", count = nrow(common_resolved)))
}

results$type <- factor(results$type, levels = c("Annual Only", "Both", "Monthly Only"))



# Create a mapping for system names to their titles
system_titles <- c(PEB = "Port Erin Bay", LZ = "Lake Zurich", NS = "North Sea")

library(ggplot2)
library(gridExtra)  # For grid.arrange
library(RColorBrewer)  # For color palettes

# Define colors from the Dark2 palette
dark2_colors <- brewer.pal(3, "Dark2")

# Map the colors to the types
my_colors <- setNames(dark2_colors, c("Annual Only", "Both", "Monthly Only"))

# Create a list to hold ggplot objects
plot_list <- list()

# Create a separate ggplot object for each system
systems <- unique(results$system)

for (s in systems) {
  system_data <- results[results$system == s, ]
  plot_list[[s]] <- ggplot(system_data, aes(x = type, y = count, fill = type)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
    scale_fill_manual(values = my_colors) +
    labs(title = system_titles[s],
         x = NULL,  # Remove the x axis label
         y = "Number of Interactions") +
    theme_minimal() +
    theme(legend.position = "none",  # Remove the legend for individual plots
          plot.title = element_text(hjust = 0.5, size = 16),  # Center and size the plot title
          axis.title.y = element_text(size = 14),  # Y axis title size
          axis.text.x = element_text(size = 15),  # X axis text size
          axis.text.y = element_text(size = 15))  # Y axis text size
}

# Use grid.arrange to plot all ggplot objects in a single plot
grid.arrange(grobs = plot_list, ncol = 1)


