
# Define a function to load data and identify resolved interactions
identify_resolved_interactions <- function(filepath) {
  load(filepath)
  resolved <- ccm_output_combined$convergence_p_val < 0.05 & ccm_output_combined$ccm_rho > ccm_output_combined$`0.95P_trials`
  ccm_output_combined[resolved, c("target", "lib")]
}
table = {}
count = 1
systems = c("Lake Zurich", "North Sea", "Port Erin Bay")
for(s in c("LZ", "NS", "PEB")){

monthly_path <- paste0(getwd(), "/Data/CCM Output/", s, "_CCM_monthly.RData")
annual_path <- paste0(getwd(), "/Data/CCM Output/", s, "_CCM_annual.RData")
load(paste0(getwd(),"/Data/Ecosystem Data/", s, "_ecosystem_data.RData"))

monthly_resolved <- identify_resolved_interactions(monthly_path)
annual_resolved <- identify_resolved_interactions(annual_path)
# Assuming taxa, monthly_resolved, and annual_resolved are already defined
# Counting occurrences in the monthly and annual data
monthly_counts <- aggregate(lib ~ target, monthly_resolved, FUN = length)
annual_counts <- aggregate(lib ~ target, annual_resolved, FUN = length)

# Renaming columns for clarity
names(monthly_counts)[2] <- "monthly_count"
names(annual_counts)[2] <- "annual_count"

# Merging the monthly and annual counts
merged_counts <- merge(monthly_counts, annual_counts, by = "target")

# Calculating the ratio of monthly to annual counts
merged_counts$ratio <- with(merged_counts, monthly_count / annual_count)

# Displaying the merged data with ratios
taxa_aggregate_df <- as.data.frame(taxa_aggregate, stringsAsFactors = FALSE)
names(taxa_aggregate_df) <- c("target", "groups_used")

# Merging the merged_counts with taxa_aggregate
final_output <- merge(merged_counts, taxa_aggregate_df, by = "target")
final_output$system <- systems[count]
final_output = final_output[,c('system', 'groups_used', 'target', 'monthly_count', 'annual_count')]
if(s == "LZ"){
  final_output = merge(final_output, taxa, by.x = "target", by.y = "id_CH")
  final_output$name <- apply(final_output[, c("order", "genus", "species", "stage")], 1, function(x) paste(na.omit(x), collapse=" "))
  final_output = final_output[,c('system', 'groups_used', 'name', 'monthly_count', 'annual_count')]
  
}
colnames(final_output) = c("System", "Functional Group", "Species", "Monthly Interactions", "Annual Interactions")
table = rbind(table,final_output)
count = count+1
}
table$`Functional Group` <- tools::toTitleCase(as.character(o$`Functional Group`))

