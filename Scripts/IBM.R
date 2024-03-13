library(parallel)
data <- list()
growth_rate = 5000
grid_size = 1000
initial_sizes <- list()
initial_sizes$nutrients = 5000
initial_sizes$plankton = 500
initial_sizes$whales = 10

eat_counts <- list()
eat_counts$nutrients = 1
eat_counts$plankton = 1
eat_counts$whales = 1000

eat_grows <- list()
eat_grows$nutrients = 2
eat_grows$plankton = 2
eat_grows$whales = eat_counts$whales * 6

eat_dies <- list()
eat_dies$nutrients = 1
eat_dies$plankton = 1
eat_dies$whales = eat_counts$whales * 3

initial_eat_count_range <- list()
initial_eat_count_range$nutrients = 1
initial_eat_count_range$plankton = 3
initial_eat_count_range$whales = 10


diet <- list()
diet$nutrients = NA
diet$plankton = "nutrients"
diet$whales = "plankton"


sizes <- list()
sizes$nutrients = 1
sizes$plankton = 30
sizes$whales = 50


speeds <- list()
speeds$nutrients = 10
speeds$plankton = 10
speeds$whales = 20

for(group in c("nutrients", "plankton", "whales")){
  curr = data.frame(replicate(3,sample(NA,initial_sizes[[group]],rep=TRUE)))
  colnames(curr) = c("x", "y", "eats")
  curr$x =runif(initial_sizes[[group]], min=1, max=grid_size)
  curr$y = runif(initial_sizes[[group]], min=1, max=grid_size)
  
  curr$eats <- lapply(1:initial_sizes[[group]], function(x) {
    sample(c(0:initial_eat_count_range[[group]]), eat_counts[[group]], replace = T)
  })
  
  data[[group]] = curr
}

simulate_eating <- function(index, radius){
    x = curr[index,1]
    y = curr[index,2]
    distances = sqrt((food[,1] - x)^2 + (food[,2] - y)^2)
    vals = which(distances < radius)
    return(vals)
}
timeseries = {}
for(i in c(1:10000)){
  print(i)
  first = T
  abundances = {}
  for(group in c("nutrients", "plankton", "whales")){
    
    curr = data[[group]]
    new_x = curr[,1] + sample(c(-speeds[[group]]:speeds[[group]]), nrow(curr),replace = T )
    new_x[new_x < 0] = 0
    new_x[new_x > grid_size] = grid_size
    new_y = curr[,2] + sample(c(-speeds[[group]]:speeds[[group]]), nrow(curr),replace = T )
    new_y[new_y < 0] = 0
    new_y[new_y > grid_size] = grid_size
    curr[,1] <- new_x
    curr[,2] <- new_y
    
   
    if(!is.na(diet[[group]])){
      food = data[[diet[[group]]]]
      new_eats = array(NA, dim = c(1, nrow(curr)))
      
      radius = sizes[[group]]
      feeds = mclapply(c(1:nrow(curr)), simulate_eating, 
                       radius = radius,
                       mc.cores = 8)
      used = {}
      for(i in c(1:length(feeds))){
        vals = feeds[[i]]
        if(sum(vals %in% used) > 0){
          vals = vals[which(-vals %in% used)]
        }
        used = c(used, vals)
        new_eats[i] = length(vals)
      }
      new_eats = as.numeric(new_eats)
      if(length(used) > 0){
        food = food[-used,]
      }
      print(paste(group, "ate", length(used), diet[[group]]))
      
      
      counts = curr$eats
      if(eat_counts[[group]] > 1){
        for(i in c(1:length(counts))){
          counts[[i]] = c(counts[[i]][2:length(counts[[i]])], as.numeric(new_eats[i]))
        }
        curr$eats = counts
        recent_consumption = unlist(lapply(counts, sum))
      }else{
        for(i in c(1:length(counts))){
          counts[[i]] = new_eats[i]
        }
        curr$eats = counts
        recent_consumption = new_eats
       }
      if(length(which(recent_consumption < eat_dies[[group]])) >0){ 
        print(paste(length(which(recent_consumption < eat_dies[[group]])), group, "starved"))
        curr = curr[-which(recent_consumption < eat_dies[[group]]),]
      }
      if(length(which(recent_consumption > eat_grows[[group]])) >0){ 
        l =  length(which(recent_consumption > eat_grows[[group]]))
        l = max(2, l)
        print(paste(l, group, "born"))
        
        curr2 = data.frame(replicate(3,sample(NA,l,rep=TRUE)))
        colnames(curr2) = c("x", "y", "eats")
        curr2$x =runif(l, min=1, max=grid_size)
        curr2$y = runif(l, min=1, max=grid_size)
        
        curr2$eats <- lapply(1:l, function(x) {
          sample(c(0:(initial_eat_count_range[[group]]) ), eat_counts[[group]], replace = T)
        })
        
        if(l == 1){
          curr = rbind(curr,  curr2[1,])
        }else{
        curr = rbind(curr,  curr2)
        }
      }
      data[[diet[[group]]]] = food
    }else{
      curr2 = data.frame(replicate(3,sample(NA,growth_rate,rep=TRUE)))
      colnames(curr2) = c("x", "y", "eats")
      curr2$x =runif(growth_rate, min=1, max=grid_size)
      curr2$y = runif(growth_rate, min=1, max=grid_size)
      
      curr2$eats <- lapply(1:growth_rate, function(x) {
        sample(c(1:10), eat_counts[[group]], replace = T)
      })
      curr=  rbind(curr, curr2)
    }
    data[[group]] = curr
    
    

  }
  
 
    
  
  
  dev.off()
  abundances = {}
  for(group in c("nutrients", "plankton", "whales")){
   abundances = c(abundances, nrow(data[[group]])) 
   curr = data[[group]]
   if(group == "nutrients"){
     plot(curr[,1], curr[,2], pch = 20, col = 'grey',
          xlim = c(0,1000), ylim = c(0,1000))
   }
   
   if(group == 'plankton'){
     points(curr[,1], curr[,2], col = rgb(.2,.6,.3), pch = 19, lwd = 3)
   }
   if(group == 'whales'){
     points(curr[,1], curr[,2], col = 'blue', pch = 19, lwd = 10)
   }
  }
  Sys.sleep(1)
  timeseries = rbind(timeseries, abundances)
  print(abundances)
  print("_______________________________")
  par(mfrow = c(3,1), mai = c(.2,.2,.2,.2) )
  if(nrow(timeseries) > 1005){
  for(i in c(1:3)){
    plot(timeseries[c(1001:nrow(timeseries)),i], type = 'l')
  }
  }
}

library(rEDM)
smooth <- function(ts, m){
  avg = array(NA, dim = c(m-1,1))
  for(sm in c(m:length(ts))){
    curr = mean(ts[(sm-m+1):sm],na.rm = T)
    avg = c(avg, curr)
  }
  return(avg)
}

vals = {}
for(sm in seq(1, 400, 20)){
  m = cbind(timeseries[,2], timeseries[,3])
  m[,2] <- smooth(m[,2], sm)
  m=m[complete.cases(m),]
  k = ccm(m, target_column = 1, lib_column = 2,
    E = 14, tp = -1, num_samples = 1, random_libs = F,
    lib_sizes = nrow(timeseries),
    lib = c(1, nrow(timeseries)),
    pred = c(1, nrow(timeseries)))
  print(k$rho)
  cs = ccf(m[,1], m[,2], plot = F)
  m = max(abs(as.numeric(cs$acf)))
  vals = c(vals, k$rho - m)
  plot(vals)
}







