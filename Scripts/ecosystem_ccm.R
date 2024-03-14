library(rEDM) #Version 0.7.5
library(parallel)
library(zoo)

# Define functions
mae <- function(x, y) {
  mean(abs(x - y), na.rm = TRUE)
}

normalize <- function(data, method = 1) {
  ts <- data[!is.na(data)]
  if (method == 1) {
    (data - min(ts)) / (max(ts) - min(ts))
  } else if (method == 2) {
    (data - mean(ts)) / sd(ts)
  }
}

smooth <- function(ts, m) {
  sapply(m:length(ts), function(sm) mean(ts[(sm - m + 1):sm], na.rm = TRUE))
}

ccm_fw <- function(target, lib, ts_matrix, interaction_matrix, num_trials = 100, E = 12, tau = 1, tp = 0, randomize = 'monthly') {
  matrix <- cbind(ts_matrix[, target], ts_matrix[, lib])
  
  # Convergence Check
  out <- ccm(matrix, target_column = 1, lib_column = 2, E = E, 
             lib_sizes = c(ceiling(nrow(matrix)/5), ceiling(nrow(matrix)/2)), 
             num_samples = 50, silent = TRUE, tau = tau, tp = tp)
  
  sample1 <- out$rho[out$lib_size == ceiling(nrow(matrix)/5)]
  sample2 <- out$rho[out$lib_size == ceiling(nrow(matrix)/2)]
  converge_significance_rho <- t.test(sample1, sample2, paired = FALSE, alternative = "less")
  
  # Main CCM Analysis
  out_main <- ccm(matrix, target_column = 1, lib_column = 2, E = E, tp = tp,
                  lib_sizes = nrow(matrix), random_libs = FALSE, 
                  num_samples = 1, silent = TRUE, tau = tau)
 
  
  # Randomization Trials
  trials_rho <- numeric(num_trials)
  for (trial in 1:num_trials) {
    shuffled_matrix <- matrix
    if (randomize == "monthly") {
      for (r in 1:12) {
        monthly_indices <- seq(r, nrow(matrix), 12)
        shuffled_matrix[monthly_indices, 2] <- shuffled_matrix[sample(monthly_indices), 2]
      }
    } else {
      shuffled_matrix[, 2] <- sample(shuffled_matrix[, 2])
    }
    
    out_trial <- ccm(shuffled_matrix, target_column = 1, lib_column = 2, E = E, tp = tp,
                     lib_sizes = nrow(matrix), random_libs = FALSE, tau = tau,
                     num_samples = 1, silent = TRUE)
    
   
    trials_rho[trial] <- out_trial$rho
  }
  
  # Preparing Output
  output <- data.frame(
    target = colnames(ts_matrix)[target],
    lib = colnames(ts_matrix)[lib],
    convergence_p_value = converge_significance_rho$p.value,
    ccm_rho = out_main$rho,
    average_trials_rho = mean(trials_rho, na.rm = TRUE),
    p95_trials_rho = quantile(trials_rho, 0.95, na.rm = TRUE),
    p99_trials_rho = quantile(trials_rho, 0.99, na.rm = TRUE)
  )
  
  return(output)
}

# Process data
process_data <- function(s, ts_matrix, interaction_matrix, n_cores, randomize, E, tau, tp, filename_suffix) {
  ccm_output <- lapply(1:ncol(ts_matrix), function(i) {
    message("Processing column: ", i)
    output <- mclapply(1:ncol(ts_matrix), ccm_fw, lib = i, ts_matrix = ts_matrix, interaction_matrix = interaction_matrix,
                       num_trials = 100, E = E, tp = tp, tau = tau, randomize = randomize, mc.cores = n_cores)
  })
  ccm_output_flat <- unlist(ccm_output, recursive = FALSE)
  
  ccm_output_combined <- do.call(rbind, ccm_output_flat)
  colnames(ccm_output_combined) <- c('target', 'lib', 'convergence_p_val', 'ccm_rho', 'average_trials', '0.95P_trials', '0.99P_trials')
  save(ccm_output_combined, file = paste0(s, "_CCM_", filename_suffix, ".RData"))
}

# Main execution
datasets <- c("LZ", "PEB", "KF", "NS")
n_cores <- detectCores()

for (s in datasets) {
  load(paste0(getwd(), "/Data/Ecosystem Data/",s,"_ecosystem_data.RData"))
  ts_matrix <- apply(ts_matrix, 2, function(x) sqrt(normalize(x, method = 1)))
  
  unique_aggregates <- unique(taxa_aggregate[, 2])
  aggregated_ts <- sapply(unique_aggregates, function(aggregate) {
    aggregate_indices <- which(taxa_aggregate[, 2] == aggregate)
    if (length(aggregate_indices) > 1) {
      rowMeans(ts_matrix[, aggregate_indices])
    } else {
      ts_matrix[, aggregate_indices]
    }
  })
  colnames(aggregated_ts) <- unique_aggregates
  aggregated_ts <- as.data.frame(aggregated_ts)
  
  if (s == 'KF') {
     process_data(s, ts_matrix, interaction_matrix, n_cores, 'random', 4, 1, -1, 'annual')
     process_data(s, aggregated_ts, interaction_matrix, n_cores, 'random', 4, 1, -1, 'aggregated')
  } else {
     process_data(s, ts_matrix, interaction_matrix, n_cores, 'monthly', 12, 1, 0, 'monthly')
     process_data(s, ts_matrix, interaction_matrix, n_cores, 'random', 4, 12, 0, 'annual')
     process_data(s, aggregated_ts, interaction_matrix, n_cores, 'monthly', 12, 1, -1, 'aggregated')
  }
}
