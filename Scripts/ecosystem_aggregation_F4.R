# Define the function for identifying significant interactions
identify_resolved_interactions <- function(filepath) {
  load(filepath)
  resolved <- ccm_output_combined$convergence_p_val < .05 & ccm_output_combined$ccm_rho > ccm_output_combined$`0.95P_trials`
  ccm_output_combined[resolved, c("target", "lib")]
}

# Main execution
datasets <- c("PEB", "LZ", "NS", "KF")
results <- list()

for (s in datasets) {
  # Load taxa_aggregate for the current dataset
  load(paste0(getwd(),"/Data/Ecosystem Data/", s, "_ecosystem_data.RData"))
  taxa_aggregate <- as.data.frame(taxa_aggregate)
  
  # Calculate the number of species in each aggregate
  species_count_in_aggregate <- table(taxa_aggregate$groups_used)
  
  # Determine the correct file path based on the dataset
  data_suffix <- if(s == "KF") "annual" else "monthly"
  filepath_non_aggregated <- paste0(getwd(),"/Data/CCM Output/", s, "_CCM_", data_suffix, ".RData")
  
  # Identify significant interactions for non-aggregated data
  significant_non_aggregated <- identify_resolved_interactions(filepath_non_aggregated)
  
  # Map the species to their aggregate groups
  significant_non_aggregated$target_group <- taxa_aggregate[match(significant_non_aggregated$target, taxa_aggregate$names), "groups_used"]
  significant_non_aggregated$lib_group <- taxa_aggregate[match(significant_non_aggregated$lib, taxa_aggregate$names), "groups_used"]
  
  # Count the occurrences of each aggregate-aggregate pair
  pair_counts <- table(significant_non_aggregated$target_group, significant_non_aggregated$lib_group)
  
  # Normalize the counts by potential interactions
  normalized_counts <- pair_counts
  for (target_group in rownames(normalized_counts)) {
    for (lib_group in colnames(normalized_counts)) {
      total_interactions <- species_count_in_aggregate[target_group] * species_count_in_aggregate[lib_group]
      normalized_counts[target_group, lib_group] <- normalized_counts[target_group, lib_group] / total_interactions
    }
  }
  
  # Store results
  results[[s]] <- normalized_counts
}


# Prepare data for plotting
plot_data <- list()

for (s in datasets) {
  # Flatten the table to a data frame
  df <- as.data.frame(as.table(results[[s]]))
  colnames(df) <- c("target_group", "lib_group", "interaction_ratio")
  
  # Load the corresponding significant aggregated data
  filepath_aggregated <- paste0(getwd(),"/Data/CCM Output/", s, "_CCM_aggregated.RData")
  significant_aggregated <- identify_resolved_interactions(filepath_aggregated)
  
  # Determine if each pair has an aggregate link
  df$has_link <- apply(df[, 1:2], 1, function(row) {
    any(significant_aggregated$target == row[1] & significant_aggregated$lib == row[2])
  })
  
  df$has_link <- ifelse(df$has_link, "Link", "No Link")
  df$System <- s  # Add a column for the system
  df = df[which(df$target_group != df$lib_group),]
  plot_data[[s]] <- df
}
# Assuming 'all_plot_data' contains the combined plot data for all datasets as previously prepared

perform_one_sided_t_tests <- function(plot_data) {
  systems <- unique(plot_data$System)
  p_values <- list()
  
  for (system in systems) {
    system_data <- subset(plot_data, System == system)
    if (length(unique(system_data$has_link)) == 2) {
      # Perform a one-sided t-test
      test_result <- wilcox.test(interaction_ratio ~ has_link, data = system_data, alternative = "greater")
      p_values[system] <- test_result$p.value
    } else {
      p_values[system] <- NA  # NA for systems where the test is not applicable
    }
  }
  
  return(p_values)
}

# Assuming all_plot_data is already prepared

# Load necessary libraries
library(ggplot2)

# Combine plot data from all datasets
all_plot_data <- do.call(rbind, plot_data)
link_counts <- table(all_plot_data$System, all_plot_data$has_link)

# Replace abbreviated names with full names in the all_plot_data dataframe
all_plot_data$System <- factor(all_plot_data$System,
                               levels = c("LZ", "PEB", "NS", "KF"),
                               labels = c("Lake Zurich", "Port Erin Bay", "North Sea", "Kelp Forest"))

# Create boxplots with individually scaled y-axes
ggplot(all_plot_data, aes(x = has_link, y = interaction_ratio, fill = has_link)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~System, scales = "free_y") +  # Individual y-axis scales for each plot
  labs(title = "", 
       x = "Functional Group Interaction", 
       y = "Fine-Scale Connectance") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 16))
t_test_results <- perform_one_sided_t_tests(all_plot_data)
link_counts
