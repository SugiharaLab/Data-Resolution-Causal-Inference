# Parameters
library(ggplot2)
library(parallel)
library(rEDM) #Version 0.7.5
n_steps <- 1000
n_predators <- 10
n_prey <- 10
output = {}
perform_trial <- function(trial) {
  
connectance <- sample(c(3:9), 1) / 10 # Connectance parameter
autocorrelation_strength <- 0.15 # Strength of negative autocorrelation

# Initial abundances (between 0 and 1)
prey_abundance <- runif(n_prey, 0, 1)
predator_abundance <- runif(n_predators, 0, 1)

# Interaction matrices
prey_interaction_matrix <- matrix(runif(n_predators * n_prey) < connectance, nrow = n_predators, ncol = n_prey) * runif(n_predators * n_prey, -1, 0)
predator_interaction_matrix <- matrix(runif(n_predators * n_prey) < connectance, nrow = n_prey, ncol = n_predators) * runif(n_predators * n_prey, 0, 1)

# Update function with negative autocorrelation
update_abundances <- function(predators, prey, prey_interaction, predator_interaction, autocorrelation_strength) {
  new_prey_abundance <- prey + rowSums(prey_interaction * predators) - autocorrelation_strength * prey
  new_predator_abundance <- predators + colSums(predator_interaction * prey) - autocorrelation_strength * predators
  
  # Apply the constraint logic
  if (any(new_prey_abundance > 1)) {
    new_prey_abundance[new_prey_abundance > 1] <- 1 / new_prey_abundance[new_prey_abundance > 1]
  }
  if (any(new_predator_abundance > 1)) {
    new_predator_abundance[new_predator_abundance > 1] <- 1 / new_predator_abundance[new_predator_abundance > 1]
  }
  
  return(list(predators = new_predator_abundance, prey = new_prey_abundance))
}

# Store results for plotting
prey_over_time <- matrix(NA, nrow = n_steps, ncol = n_prey)
predators_over_time <- matrix(NA, nrow = n_steps, ncol = n_predators)

# Simulation loop
for (t in 1:n_steps) {
  result <- update_abundances(predator_abundance, prey_abundance, prey_interaction_matrix, predator_interaction_matrix, autocorrelation_strength)
  predator_abundance <- result$predators
  prey_abundance <- result$prey
  
  # Store the results
  prey_over_time[t, ] <- prey_abundance
  predators_over_time[t, ] <- predator_abundance
}

m = cbind(rowSums(prey_over_time), rowSums(predators_over_time))[200:length(rowSums(predators_over_time)),]

for(i in c(1:ncol(m))){
  m[,i] = m[,i] * sample(c(900:1100), nrow(m), replace=  T)/1000
}

val = ccm(m, target_column = 1, lib_column = 2, random_libs = F,
    num_samples = 1, lib_sizes = nrow(m), E = 10, tp = -1, exclusion_radius = 100,
    )$rho
return(c(connectance, val))
}

# Run the trials in parallel
output <- mclapply(1:500, perform_trial, mc.cores = detectCores() - 1)

# Convert the list to a matrix
output_matrix <- do.call(rbind, output)

# Convert the matrix to a data frame
output_df <- as.data.frame(output_matrix)
names(output_df) <- c("Group", "Value")

# Convert the first column to a factor
output_df$Group <- as.factor(output_df$Group)
# Load the ggplot2 package

# Convert the first column to a factor with a more descriptive name
# output_df$Group <- factor(output_df$Group, labels = paste("C =", unique(output_df$Group)))

# Create the boxplot using ggplot2
ggplot(output_df, aes(x = Group, y = Value)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  labs(title = "",
       x = "Connectance",
       y = "Resolved Aggregate Interaction Strength") +
  theme_minimal() +  # Minimal theme for a clean look
  theme(axis.title = element_text(size = 12, face = "bold"),  # Customize axis labels
        axis.text = element_text(size = 10),  # Customize axis text
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))  # Center and style the title
