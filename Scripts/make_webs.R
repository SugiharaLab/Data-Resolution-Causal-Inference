# Load required libraries
library(plotrix)
library(RColorBrewer)
library(shape)
library(stringr)
library(corrplot)

# Set graphical parameters
par(mfrow = c(2, 2), mai = c(1, 1, 1, 1) * 0.65, cex = 1.75)

# Titles for the plots
titles = rbind(c("NS", "North Sea"),
               c("PEB", "Port Erin Bay"),
               c("LZ", "Lake Zurich"),
               c("KF", "Kelp Forest"))

# Function to draw arrows with offset for bidirectional arrows
draw_arrow <- function(x1, y1, x2, y2, radius, width = 0.1, length = 0.1, lwd = 1, col = 'black', offset_factor = 0.1) {
  if (x2 != x1 || y2 != y1) {
    # Calculate direction vector
    dir_vector <- c(x2 - x1, y2 - y1)
    dir_vector <- dir_vector / sqrt(sum(dir_vector^2))
    
    # Calculate intersection points with the circle
    intersect_start <- c(x1, y1) + radius * dir_vector
    intersect_end <- c(x2, y2) - radius * dir_vector
    
    # Apply additional offset to avoid overlap in bidirectional arrows
    offset = offset_factor * radius
    intersect_start <- intersect_start + offset * dir_vector
    intersect_end <- intersect_end - offset * dir_vector
    
    # Draw arrow
    Arrows(intersect_start[1], intersect_start[2], intersect_end[1], intersect_end[2], 
           lwd = lwd, col = col, arr.type = "triangle", arr.width = width, arr.length = length)
  }
}

make_web <- function(taxa_aggregate, arrow_details, r = .2, names = F, dots = T){
  aggregate_count = length(unique(taxa_aggregate[,2]))
  colors = brewer.pal(aggregate_count, 'Set3')
  aggregate_coordinates = {}
  for(i in c(0: (aggregate_count-1) )){
    aggregate_coordinates = rbind(aggregate_coordinates, 
                                  c(cos(2*pi*i / aggregate_count ),sin(2*pi*i / aggregate_count ) ))
  }
  rownames(aggregate_coordinates)= as.character(unique(taxa_aggregate[,2]))
  species_names = {}
  species_coordinates = {}
  for(aggregate in unique(taxa_aggregate[,2])){
    curr = taxa_aggregate[taxa_aggregate[,2] == aggregate,]
    if( !is.null(nrow(curr))){
      species_names = c(species_names, as.character(curr[,1]))
      species_count = nrow(curr)
      curr_agg_coords = aggregate_coordinates[as.character(curr[1,2]),]
      for(i in c(0 : (species_count-1) )){
        species_coordinates = rbind(species_coordinates, 
                                    c(cos(2*pi*i / species_count ) * r + curr_agg_coords[1],sin(2*pi*i / species_count ) * r + curr_agg_coords[2]))
      }
    }else{
      curr = taxa_aggregate[taxa_aggregate[,2] == aggregate,]
      
      species_names = c(species_names, as.character(curr[1]))
      species_count = 1
      curr_agg_coords = aggregate_coordinates[as.character(curr[2]),]
      species_coordinates = rbind(species_coordinates, 
                                  c( curr_agg_coords[1],curr_agg_coords[2]))
      
    }
    
    
  }
  rownames(species_coordinates) = species_names
  t = titles[titles[,1]==s,2]
  if(dots){
    plot(species_coordinates[,1], species_coordinates[,2],bty="n", yaxt = 'n', xaxt = 'n',
         xlim = c(-1,1)*1.5,ylim = c(-1,1)*1.5, xlab = '', ylab = '', pch = 19, main = t, 
         col = 'white', asp = 1)
  }else{
    plot(0, 0, col = 'white', bty="n", yaxt = 'n', xaxt = 'n',asp = 1,
         xlim = c(-1,1)*1.5,ylim = c(-1,1)*1.5, xlab = '', ylab = '', main = t)
  }
  for (i in c(1:nrow(aggregate_coordinates))) {
    cols = c(col2rgb(colors[i]) / 255, 0.5)
    draw.circle(aggregate_coordinates[i, 1], aggregate_coordinates[i, 2], r * 1.5,
                col = rgb(cols[1], cols[2], cols[3], cols[4])
    )
    
    if (names) {
      n = rownames(aggregate_coordinates)[i]
      n = gsub(" ", "\n", n)
      # Calculate radial position for text based on the length of the string
      string_length_factor = 1.1 + 0.03 * nchar(n) # Adjust the factor based on string length
      radial_offset = r * 1.5 * string_length_factor
      text_x = aggregate_coordinates[i, 1] * (1 + radial_offset)
      text_y = aggregate_coordinates[i, 2] * (1 + radial_offset)
      text(text_x, text_y, str_to_title(n), font = 2, cex = 0.8)
    }
  }
  if(dots){
    points(species_coordinates[,1], species_coordinates[,2],bty="n", yaxt = 'n', xaxt = 'n',
           xlim = c(-1,1)*1.5,ylim = c(-1,1)*1.5, xlab = '', ylab = '', pch = 19)
  }
  coordinates = rbind(aggregate_coordinates, species_coordinates)
  if( !is.null(nrow(arrow_details) )){
    for(i in c(1:nrow(arrow_details))){
      if(arrow_details[i,5] == "aggregate"){
        draw_arrow(as.numeric(coordinates[arrow_details[i,1],1]),
                   as.numeric(coordinates[arrow_details[i,1],2]),
                   as.numeric(coordinates[arrow_details[i,2],1]),
                   as.numeric(coordinates[arrow_details[i,2],2]),
                   r * 2, lwd = as.numeric(arrow_details[i,4]), 
                   col = arrow_details[i,3] )
      }else{
        draw_arrow(as.numeric(coordinates[arrow_details[i,1],1]),
                   as.numeric(coordinates[arrow_details[i,1],2]),
                   as.numeric(coordinates[arrow_details[i,2],1]),
                   as.numeric(coordinates[arrow_details[i,2],2]),
                   r *0.15, lwd = as.numeric(arrow_details[i,4]), 
                   col = arrow_details[i,3],
                   width = .07, length = .1)
      }
    }
  }
  
  
  
  
  
  
}


for(s in c("PEB","LZ", "NS", "KF")){
  filepath_aggregated <- paste0(getwd(),"/Data/CCM Output/", s, "_CCM_aggregated.RData")
  load(paste0(getwd(),"/Data/Ecosystem Data/", s, "_ecosystem_data.RData"))
  load(filepath_aggregated)
  arrow_details = {} 
  #FOOD WEB
  for(i in c(1:nrow(directed_fw))){
    for(j in c(1:nrow(directed_fw))){
      if(directed_fw[i,j] == 1){
        arrow_details = rbind(arrow_details, c(colnames(directed_fw)[i], colnames(directed_fw)[j], 'red', 10 ,"aggregate"))
      }
    }
  }

  #AGGREGATE CAUSAL
  for(i in c(1:nrow(ccm_output_combined))){
    if(as.numeric(ccm_output_combined[i,]$convergence_p_val) < 0.05 &
       as.numeric(ccm_output_combined[i,]$ccm_rho) >  as.numeric(ccm_output_combined[i,][['0.95P_trials']]) ){
        arrow_details = rbind(arrow_details, c(ccm_output_combined$target[i], ccm_output_combined$lib[i], 'blue', 2.5, "aggregate"))

    }
  }



  # data_suffix <- if(s == "KF") "annual" else "monthly"
  # filepath_non_aggregated <- paste0(getwd(), "/Data/CCM Output/", s, "_CCM_", data_suffix, ".RData")
  # load(filepath_non_aggregated)
  # #HIGH RES CAUSAL
  # for(i in c(1:nrow(ccm_output_combined))){
  #   if(as.numeric(ccm_output_combined[i,]$convergence_p_val) < 0.05 &
  #      as.numeric(ccm_output_combined[i,]$ccm_rho) >  as.numeric(ccm_output_combined[i,][['0.95P_trials']]) ){
  #   
  #     arrow_details = rbind(c(ccm_output_combined$target[i], ccm_output_combined$lib[i], rgb(.2,.2,.2,.5), 1, "species"),arrow_details)
  #   }
  # }
  # arrow_details = unique(arrow_details)
  # arrow_details = arrow_details[arrow_details[,1] != arrow_details[,2],]
  colnames(arrow_details) = c("from", "to","color", "lwd", "type")
  
  
 
  make_web(taxa_aggregate, arrow_details, names = T, dots = T)
 
}





