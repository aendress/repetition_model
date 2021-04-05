rm (list = ls())

library (tidyr)
library(plyr)
library (ggplot2)
library (ggthemes)

source ("/Users/endress/R.ansgar/tt.R")

# Set parameters

# In the sequential case, this is the number of neurons per layer, in the simultaneous case one of the 
# dimensions of the layers
N.FEATURES <- 10
# ignored for the sequential case
N.LOCATIONS <- 5

NOISE.ACTIVATION.SPREAD <- 0
NOISE.ACTIVATION.SPREAD.SEQ <- seq (0, .25, len=50)
ADD.NOISE.BEYOND.INPUT <- TRUE

N.SUBJ <- 50

STANDARD.WEIGHT <- 1
NOISE.WEIGHT.SPREAD <- 0

THRESHOLD.OUTPUT <- .5
TONIC.INHIBITION <- 1


# Multiple plot function
# From http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# 
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

add.noise <- function (activation.matrix,
                       noise.activation.spread = NOISE.ACTIVATION.SPREAD,
                       noise.activtion.m = 0){
  
  if (is.matrix (activation.matrix)){
    activation.matrix.dim <- dim (activation.matrix)
  } else {
    if (is.numeric (activation.matrix) && length (activation.matrix) == 2){
      activation.matrix.dim <- activation.matrix  
    } else {
      if (!is.numeric (activation.matrix) || !is.vector(activation.matrix))
        stop ("First argument needs to be a matrix, the dimensions of a matrix, or a vector")
    }
  }
  
  if (is.matrix(activation.matrix)){
    # Input is a matrix
    return (activation.matrix +     
              matrix (generate.random.numbers.around.target (prod (activation.matrix.dim),
                                                             noise.activtion.m, noise.activation.spread),
                      nrow = activation.matrix.dim[1],
                      ncol = activation.matrix.dim[2]))
  } else {
    if (length (activation.matrix) == 2){
      # input are dimensions of a matrix
      return (matrix (generate.random.numbers.around.target (prod (activation.matrix.dim),
                                                             noise.activtion.m, noise.activation.spread),
                      nrow = activation.matrix.dim[1],
                      ncol = activation.matrix.dim[2]))
    } else {
      # input is a vector
      return (activation.matrix +
                generate.random.numbers.around.target (length (activation.matrix),
                                                       noise.activtion.m, noise.activation.spread))
    }
  }
}  

generate.random.numbers.around.target <- function (n, target.val, spread,
                                                   noise.distribution = c("normal", "uniform")){
  # Choose lower and upper bound so that the target value is as centered as possible,
  # but the range doesn't fall out of the interval [0, 1]
  
  noise.distribution = match.arg(noise.distribution)
  
  if (spread == 0)
    return (rep(target.val, n))
  
  if (noise.distribution == "normal"){
    return (rnorm (n, target.val, spread))
  } else {
    # uniform distribution
    value.range <- target.val + c(-spread, spread) / 2 
    
    if (target.val > .5){
      # Shift leftwards
      if (value.range[2] > 1)
        value.range <- value.range - (value.range[2] - 1)
    }
    
    if (target.val < .5){
      # Shift rightwards
      if (value.range[1] < 0)
        value.range <- value.range - value.range[1]
    }
    
    if (value.range[1] < 0)
      stop ("Lower bound of random values can cross 0, please reduce spread.")
    
    if (value.range[2] > 1)
      stop ("Upper bound of random values can cross 1, please reduce spread.")
    
    message("Using uniform distribution between ",  value.range[1], " and ", value.range[2], " for target value ", target.val, "." )
    runif (n, value.range[1], value.range[2])
  }
}

restrict.value.range <- function (m, min.val = 0, max.val = 1){
  
  m.restricted <- m
  
  m.restricted[m.restricted<min.val] <- min.val
  m.restricted[m.restricted>max.val] <- max.val
  
  m.restricted  
  
}

detect.simulation.type <- function (result.list,
                                    verbose = TRUE){
  
  if (all (unlist (lapply (result.list[grep ("activation", names (result.list))], 
                           is.numeric))) &&
      all (unlist(lapply (result.list[grep ("activation", names (result.list))], 
                          function (X) is.null(dim(X)))))){
    
    if (verbose)
      message ("Sequential simulation detected")
    
    return ("sequential")
  }
  
  if (all (unlist (lapply (result.list[grep ("activation", names (result.list))], 
                           is.matrix))) &&
      all (unlist(lapply (result.list[grep ("activation", names (result.list))], 
                          function (X) length(dim(X)) == 2)))){
    # This is the simultaneous case
    if (verbose)
      message ("Simultaneous simulation detected")
    
    return ("simultaneous")
  }
  
  stop ("We could not identify the simulation type.")
}

detect.input.list.type <- function (input.list,
                                    verbose = TRUE){
  if (is.list (input.list)){
    if (verbose){
      message ("Input list detected")
    }
    
    if (is.matrix (input.list[[length(input.list)]])){
      # This is a noiseless input matrix
      return ("noiseless.matrix")
    } else {
      stop ("input.list must be a list of matrices or the coordinates of the target.")
    }
  }
  
  if (is.numeric(input.list)){
    # These are the directly specified target coordinates
    if (verbose){
      message ("Input coordinates detected")
    }
    return ("input.coordinates")
  }
  
  stop ("We could not identify the input.list type.")
}

extract.target.activation.from.copy.layer <- function (result.list,
                                                       input.list,
                                                       return.type = c("Z", "raw"),
                                                       print.mean = FALSE,
                                                       verbose = TRUE){
  
  # Return raw score or Z value (default)
  return.type <- match.arg (return.type)
  
  simulation.type <- detect.simulation.type (result.list,
                                             verbose)
  
  if (simulation.type == "sequential"){
    
    if (return.type == "raw"){
      return (result.list$activation.copy.layer[which (input.list[[length(input.list)]]> 0)])
    } else {
      
      if (sd(result.list$activation.copy.layer) == 0)
        return (0)
      
      return (scale(result.list$activation.copy.layer)[which (input.list[[length(input.list)]]> 0)])
    }
    
  } else {
    # Simultaneous case
    
    if (print.mean)
      message (mean(result.list$activation.copy.layer))
    
    activation.copy.layer.scaled <- result.list$activation.copy.layer - 
      mean(result.list$activation.copy.layer)
    
    activation.copy.layer.sd <- sd(activation.copy.layer.scaled)
    if (activation.copy.layer.sd == 0)
      return (0)
    activation.copy.layer.scaled <- activation.copy.layer.scaled / activation.copy.layer.sd 
    
    input.list.type <- detect.input.list.type (input.list,
                                               verbose)
    if (input.list.type == "noiseless.matrix"){
      target.input <- input.list[[length(input.list)]]     
      
      if (return.type == "raw"){
        return (result.list$activation.copy.layer[which (target.input > 0, arr.ind = TRUE)])
      } else {
        return (activation.copy.layer.scaled[which (target.input > 0, arr.ind = TRUE)])
      }
    } else {
      # These are the directly specified target coordinates
      target.input <- input.list
      if (return.type == "raw"){
        return (result.list$activation.copy.layer[target.input[1], target.input[2]])
      } else {
        return (activation.copy.layer.scaled[target.input[1], target.input[2]])
      }
    }
  }
}

extract.target.activation.from.output.layer <- function (result.list,
                                                         input.list,
                                                         print.mean = FALSE,
                                                         verbose = TRUE){
  
  
  simulation.type <- detect.simulation.type (result.list,
                                             verbose)
  
  if (simulation.type == "sequential"){
    
    return (result.list$activation.output.layer[which (input.list[[length(input.list)]]> 0)])
    
  } else {
    # Simultaneous case
    
    input.list.type <- detect.input.list.type (input.list,
                                               verbose)
    if (input.list.type == "noiseless.matrix"){
      target.input <- input.list[[length(input.list)]]     
      return (result.list$activation.output.layer[which (target.input> 0, arr.ind = TRUE)])
    } else {
      # These are the directly specified target coordinates
      target.input <- input.list
      return (result.list$activation.output.layer[target.input[1], target.input[2]])
    }
  }
}


make.plot <- function (results.long,
                       measure.var = "activation",
                       group.vars = c("noise.spread", "sequence_type"),
                       facet.var = "sequence_type",
                       x.var = "noise.spread",
                       xlab = x.var,
                       ylab  = measure.var,
                       legend.position = "bottom",
                       facet.grid.vars = NULL,
                       show.data.points = FALSE){
  
  results.long.m <- summarySE(results.long, 
                              measurevar=measure.var, 
                              groupvars=c("noise.spread", "sequence_type", facet.grid.vars))
  
  results.long.m$sequence_type <- mapvalues(results.long.m$sequence_type, 
                                            from = levels(results.long.m$sequence_type), 
                                            to = sub ("^.*?\\.", "", 
                                                      levels(results.long.m$sequence_type), 
                                                      perl = TRUE))
  results.long.m$sequence_type <- mapvalues(results.long.m$sequence_type, 
                                            from = levels(results.long.m$sequence_type), 
                                            to = gsub ("\\.", " ", 
                                                       levels(results.long.m$sequence_type), 
                                                       perl = TRUE))
  
  
  my.plot <- ggplot(results.long.m, 
                    aes_string(x = x.var, 
                               y = measure.var, 
                               colour = facet.var))
  
  my.plot <- my.plot + 
    theme_minimal()
  #  theme_igray()
  #  theme_few()  
  # theme_fivethirtyeight()
  
  if (measure.var == "activation"){
    my.plot <- my.plot +
      geom_ribbon(aes(ymin = activation-se,
                      ymax = activation+se), 
                  alpha = 0.2)
    #geom_errorbar(aes(ymin=activation-se, ymax=activation+se, group="sequence_type"), width=.02) +
  }
  if (measure.var == "correct"){
    my.plot <- my.plot +
      geom_ribbon(aes(ymin = correct-se,
                      ymax = correct+se), 
                  alpha = 0.2)
    #geom_errorbar(aes(ymin=activation-se, ymax=activation+se, group="sequence_type"), width=.02) +
  }
  
  my.plot <- my.plot +
    geom_line() 
  
  if (show.data.points){
    my.plot <- my.plot +
      geom_point()
  }
  
  my.plot <- my.plot +
    xlab(xlab) 
  
  if (is.null((facet.grid.vars))){
    my.plot <- my.plot +
      ylab (ylab)
  } else {
    my.plot <- my.plot +
      theme(axis.title.y = element_blank())
  }
  
  if (!is.null (facet.grid.vars)){
    # my.plot <-  my.plot +
    #   facet_grid(reformulate(facet.grid.vars[1],facet.grid.vars[2]))
    
    if (length (facet.grid.vars) == 2){
      my.plot <-  my.plot +
        facet_grid(reformulate(facet.grid.vars[1],facet.grid.vars[2]))
    } else {
      my.plot <-  my.plot +
        facet_grid(reformulate(facet.grid.vars[1]))
    }
  }
  
  my.plot <-  my.plot +
    theme(legend.position = legend.position,
          legend.title = element_blank(),
          legend.text = element_text(size=11),
          #legend.position='none',
          legend.background = element_rect(fill = "transparent")) +
    guides(col = guide_legend(ncol = 2, byrow=TRUE))
  
  my.plot <-  my.plot +
    theme(strip.text = element_text(size = my.plot$theme$legend.text$size + 2), 
          axis.title.x = element_text(size = my.plot$theme$legend.text$size + 2), 
          axis.text = element_text(size = my.plot$theme$legend.text$size))
  
  my.plot <- my.plot + 
    theme (plot.margin = margin(.25, 1, .25, 1, "lines"))
  
  # my.plot <- my.plot +     
  #   theme (panel.background = element_rect(fill = "white"),
  #          plot.background = element_rect(
  #            fill = "grey90",
  #            colour = "black",
  #            size = 1))
  # 
  
  my.plot
}


# Sequential repetition detection

run.simulation.sequential <- function (input.list,
                                       n.features = N.FEATURES,
                                       noise.activation.spread = NOISE.ACTIVATION.SPREAD,
                                       standard.weight = STANDARD.WEIGHT, 
                                       noise.weight.spread = NOISE.WEIGHT.SPREAD, 
                                       threshold.output = THRESHOLD.OUTPUT, 
                                       tonic.inhibition = TONIC.INHIBITION){
  
  # Set up layers
  
  # receives input
  activation.input.layer <- generate.random.numbers.around.target (n.features, 
                                                                   0, noise.activation.spread)
  
  # linear neurons receiving input from input layer and inhibition layer
  activation.copy.layer <- generate.random.numbers.around.target (n.features, 
                                                                  0, noise.activation.spread)
  weights.copy.from.input <- generate.random.numbers.around.target (n.features, 
                                                                    standard.weight, noise.weight.spread)
  weights.copy.from.inhibition <- -generate.random.numbers.around.target (n.features, 
                                                                          standard.weight, noise.weight.spread)
  
  # linear neurons receiving input from input layer
  activation.inhibition.layer <- generate.random.numbers.around.target (n.features, 
                                                                        tonic.inhibition, noise.activation.spread)
  weights.inhibition.from.input <- -generate.random.numbers.around.target (n.features, 
                                                                           standard.weight, noise.weight.spread)
  
  
  # thresholded neurons receiving input from copy layer
  activation.output.layer <- generate.random.numbers.around.target (n.features, 
                                                                    0, noise.activation.spread)
  weights.output.from.copy <- generate.random.numbers.around.target (n.features, 
                                                                     standard.weight, noise.weight.spread)
  
  
  outputs <- c()
  for (current.input in input.list){
    
    activation.input.layer <- current.input
    activation.input.layer <- add.noise(activation.input.layer,
                                        noise.activation.spread = noise.activation.spread)
    
    activation.copy.layer <- (weights.copy.from.input * activation.input.layer) + 
      (weights.copy.from.inhibition * activation.inhibition.layer)
    if (ADD.NOISE.BEYOND.INPUT)
      activation.copy.layer <- add.noise (activation.copy.layer,
                                          noise.activation.spread = noise.activation.spread)
    activation.copy.layer <- restrict.value.range (activation.copy.layer)
    
    if (ADD.NOISE.BEYOND.INPUT){
      activation.inhibition.layer <- generate.random.numbers.around.target (n.features, 
                                                                            tonic.inhibition, noise.activation.spread)
    } else {
      activation.inhibition.layer <- generate.random.numbers.around.target (n.features, 
                                                                            tonic.inhibition, 0)
    }
    activation.inhibition.layer <- activation.inhibition.layer +  
      weights.inhibition.from.input * activation.input.layer
    activation.inhibition.layer <- restrict.value.range (activation.inhibition.layer)
    
    activation.output.layer <- 1 * ((weights.output.from.copy * activation.copy.layer) >= threshold.output)
    outputs <- cbind(outputs, activation.output.layer)
    
  }
  
  colnames (outputs) <- paste ("item", 1:length(input.list), sep="")
  
  list (activation.output.layer = activation.output.layer,
        output = outputs,
        activation.copy.layer = activation.copy.layer,
        activation.inhibition.layer = activation.inhibition.layer)
}


#############################################


# Simultaneous presentation

# # This version of the function used just one inhibitory population per feature. This doesn't work, as the input
# # from the activation layer is either diluted so that inhibition is hardly reduced (e.g., if the input to the
# # inhibtion layer is the average activation of a feature across locations from the input layer) or sufficiently
# # strong that a single occurrence in the activation layer can completely shut down inhibition.
# # Disinhibition thus must be spatially specific.
# run.simulation.simultaneous.not.topographic.inhibition <- function (input.matrix,
#                                          n.features = N.FEATURES,
#                                          n.locations = N.LOCATIONS,
#                                          noise.activation.spread = NOISE.ACTIVATION.SPREAD,
#                                          standard.weight = STANDARD.WEIGHT,
#                                          noise.weight.spread = NOISE.WEIGHT.SPREAD,
#                                          threshold.output = THRESHOLD.OUTPUT,
#                                          tonic.inhibition = TONIC.INHIBITION){
# 
#   # Setup layers
# 
#   # receives input
#   activation.input.layer <- matrix (generate.random.numbers.around.target (n.features * n.locations,
#                                            0, noise.activation.spread),
#                                     nrow = n.features,
#                                     ncol = n.locations)
# 
#   # linear neurons receiving input from input layer and inhibition layer
#   activation.copy.layer <- matrix (generate.random.numbers.around.target (n.features * n.locations,
#                                           0, noise.activation.spread),
#                                    nrow = n.features,
#                                    ncol = n.locations)
# 
#   weights.copy.from.input <- matrix (generate.random.numbers.around.target (n.features * n.locations,
#                                             standard.weight, noise.weight.spread),
#                                      nrow = n.features,
#                                      ncol = n.locations)
# 
#   weights.copy.from.inhibition <- -generate.random.numbers.around.target (n.features,
#                                           standard.weight, noise.weight.spread)
# 
#   # linear neurons receiving input from input layer
#   # Inhibition is NOT topographically organized, and just distinguishes features but not locations
#   activation.inhibition.layer <- generate.random.numbers.around.target (n.features,
#                                         tonic.inhibition, noise.activation.spread)
#   weights.inhibition.from.input <- -generate.random.numbers.around.target (n.features,
#                                           standard.weight, noise.weight.spread)
# 
# 
#   # thresholded neurons receiving input from copy layer
#   activation.output.layer <- matrix (generate.random.numbers.around.target (n.features * n.locations,
#                                             0, noise.activation.spread),
#                                      nrow = n.features,
#                                      ncol = n.locations)
#   weights.output.from.copy <- matrix (generate.random.numbers.around.target (n.features * n.locations,
#                                              standard.weight, noise.weight.spread),
#                                       nrow = n.features,
#                                       ncol = n.locations)
# 
# 
#   # Now run the simulation
#   activation.input.layer <- input.matrix
# 
#   activation.inhibition.layer <- generate.random.numbers.around.target (n.features,
#                                         tonic.inhibition, noise.activation.spread)
#   # Locations are in columns, so average over locations/columns
#   activation.inhibition.layer <- activation.inhibition.layer +
#     weights.inhibition.from.input * apply (activation.input.layer, 1, mean)
#   activation.inhibition.layer <- activation.inhibition.layer * (activation.inhibition.layer > 0)
# 
#   activation.copy.layer <- (weights.copy.from.input * activation.input.layer)
#   # Features are in rows, so apply sweep row-wise
#   activation.copy.layer <- sweep (activation.copy.layer,
#                                   1,
#                                   weights.copy.from.inhibition * activation.inhibition.layer,
#                                   "+")
#   activation.copy.layer <- activation.copy.layer * (activation.copy.layer > 0)
# 
#   #activation.output.layer <- 1 * ((weights.output.from.copy * activation.copy.layer) >= threshold.output)
#   activation.output.layer <- 1 * ((weights.output.from.copy * activation.copy.layer) )
#   activation.output.layer
# }


run.simulation.simultaneous.topographic.inhibition <- function (input.matrix.list,
                                                                n.features = N.FEATURES,
                                                                n.locations = N.LOCATIONS,
                                                                noise.activation.spread = NOISE.ACTIVATION.SPREAD,
                                                                standard.weight = STANDARD.WEIGHT,
                                                                noise.weight.spread = NOISE.WEIGHT.SPREAD,
                                                                threshold.output = THRESHOLD.OUTPUT,
                                                                tonic.inhibition = TONIC.INHIBITION){
  
  # Setup layers
  
  # receives input
  activation.input.layer <- add.noise (c(n.features,n.locations),
                                       noise.activation.spread = noise.activation.spread)
  
  # linear neurons receiving input from input layer and inhibition layer
  if (ADD.NOISE.BEYOND.INPUT){
    activation.copy.layer <- matrix (generate.random.numbers.around.target (n.features * n.locations,
                                                                            0, noise.activation.spread),
                                     nrow = n.features,
                                     ncol = n.locations)
  } else {
    activation.copy.layer <- matrix (generate.random.numbers.around.target (n.features * n.locations,
                                                                            0, 0),
                                     nrow = n.features,
                                     ncol = n.locations)
    
  }
  
  weights.copy.from.input <- matrix (generate.random.numbers.around.target (n.features * n.locations,
                                                                            standard.weight, noise.weight.spread),
                                     nrow = n.features,
                                     ncol = n.locations)
  
  weights.copy.from.inhibition <- -matrix (generate.random.numbers.around.target (n.features * n.locations,
                                                                                  standard.weight, noise.weight.spread),
                                           nrow = n.features,
                                           ncol = n.locations)
  
  # linear neurons receiving input from input layer
  # Inhibition is also topographically organized, and distinguishes both features but and locations
  # as the commented out model about doesn't work for the reason 
  # Each feature in the input layer disinhibits that feature across locations, and inhibits the other 
  # features, again across locations
  activation.inhibition.layer <- matrix (generate.random.numbers.around.target (n.features * n.locations,
                                                                                tonic.inhibition, noise.activation.spread),
                                         nrow = n.features,
                                         ncol = n.locations)
  weights.inhibition.from.input <- array (generate.random.numbers.around.target (n.features * n.locations * n.features * n.locations,
                                                                                 standard.weight, noise.weight.spread),
                                          dim=c(n.features, n.locations, n.features, n.locations),
                                          dimnames=list(paste("input.feature", 1:n.features, sep="."),
                                                        paste("input.location", 1:n.locations, sep="."),
                                                        paste("inhibition.feature", 1:n.features, sep="."),
                                                        paste("inhibition.location", 1:n.locations, sep="."))) 
  for (i.feat in 1:n.features){
    weights.inhibition.from.input[i.feat,,i.feat,] <- -weights.inhibition.from.input[i.feat,,i.feat,]    
  }  
  
  
  # thresholded neurons receiving input from copy layer
  activation.output.layer <- matrix (generate.random.numbers.around.target (n.features * n.locations,
                                                                            0, noise.activation.spread),
                                     nrow = n.features,
                                     ncol = n.locations)
  weights.output.from.copy <- matrix (generate.random.numbers.around.target (n.features * n.locations,
                                                                             standard.weight, noise.weight.spread),
                                      nrow = n.features,
                                      ncol = n.locations)
  
  
  # Now run the simulation
  
  # If the input is just a single matrix (i.e., the simultaneous case, 
  # put it into a list)
  if (!is.list (input.matrix.list)){
    if (is.matrix(input.matrix.list))
      input.matrix.list <- list (input.matrix.list)
    else
      stop ("First argument must be a matrix or a list of matrices")
  }
  
  for (current.input.matrix in input.matrix.list){
    activation.input.layer <- current.input.matrix
    activation.input.layer <- add.noise (activation.input.layer,
                                         noise.activation.spread = noise.activation.spread)
    
    if (ADD.NOISE.BEYOND.INPUT){
      activation.inhibition.layer <- matrix (generate.random.numbers.around.target (n.features * n.locations,
                                                                                    tonic.inhibition, noise.activation.spread),
                                             nrow = n.features,
                                             ncol = n.locations)
    } else {
      activation.inhibition.layer <- matrix (generate.random.numbers.around.target (n.features * n.locations,
                                                                                    tonic.inhibition, 0),
                                             nrow = n.features,
                                             ncol = n.locations)
    }
    
    for (i.feat in 1:n.features){
      for (i.loc in 1:n.locations){
        activation.inhibition.layer <- activation.inhibition.layer +
          weights.inhibition.from.input[i.feat,i.loc,,] * activation.input.layer[i.feat,i.loc]
      }
    }
    activation.inhibition.layer <- restrict.value.range (activation.inhibition.layer)
    
    activation.copy.layer <- (weights.copy.from.input * activation.input.layer) + 
      weights.copy.from.inhibition * activation.inhibition.layer
    
    activation.copy.layer <- restrict.value.range (activation.copy.layer)
    
    activation.output.layer <- 1 * ((weights.output.from.copy * activation.copy.layer) >= threshold.output)
  }
  
  
  list (activation.output.layer =  activation.output.layer,
        activation.copy.layer = activation.copy.layer,
        activation.inhibition.layer = activation.inhibition.layer)
  
}



run.simulation.combined <- function (input.matrix.list,
                                     n.features = N.FEATURES,
                                     n.locations = N.LOCATIONS,
                                     noise.activation.spread = NOISE.ACTIVATION.SPREAD,
                                     standard.weight = STANDARD.WEIGHT,
                                     noise.weight.spread = NOISE.WEIGHT.SPREAD,
                                     threshold.output = THRESHOLD.OUTPUT,
                                     tonic.inhibition = TONIC.INHIBITION){
  
  # Setup layers
  
  # receives input
  activation.input.layer <- add.noise(c(n.features, n.locations),
                                      noise.activation.spread = noise.activation.spread)
  
  # linear neurons receiving input from input layer and inhibition layer
  if (ADD.NOISE.BEYOND.INPUT){
    activation.copy.layer <- add.noise(c(n.features, n.locations),
                                       noise.activation.spread = noise.activation.spread)
  } else {
    activation.copy.layer <- add.noise(c(n.features, n.locations),
                                       noise.activation.spread = 0)
  }
  weights.copy.from.input <- matrix (generate.random.numbers.around.target (n.features * n.locations,
                                                                            standard.weight, noise.weight.spread),
                                     nrow = n.features,
                                     ncol = n.locations)
  
  weights.copy.from.inhibition <- -matrix (generate.random.numbers.around.target (n.features * n.locations,
                                                                                  standard.weight, noise.weight.spread),
                                           nrow = n.features,
                                           ncol = n.locations)
  
  # Self-disinhibition
  # This layer receive excitatory 1:1 input from the input layer
  # and provides inhibitory 1:1 input to the inhibition layer. It acts
  # as a delay element for self-disinhibition
  if (ADD.NOISE.BEYOND.INPUT)
    activation.selfdisinhibition.layer <- add.noise(c(n.features, n.locations),
                                                    noise.activation.spread = noise.activation.spread)
  else 
    activation.selfdisinhibition.layer <- add.noise(c(n.features, n.locations),
                                                    noise.activation.spread = 0)
  
  weights.selfdisinhibition.from.input <- matrix (generate.random.numbers.around.target (n.features * n.locations,
                                                                                         standard.weight, noise.weight.spread),
                                                  nrow = n.features,
                                                  ncol = n.locations)
  
  # linear neurons receiving input from input layer
  # Inhibition is also topographically organized, and distinguishes both features but and locations
  # as the commented out model about doesn't work for the reason 
  # Each feature in the input layer disinhibits that feature across locations, and inhibits the other 
  # features, again across locations
  if (ADD.NOISE.BEYOND.INPUT)
    activation.inhibition.layer <- add.noise(c(n.features, n.locations),
                                             noise.activtion.m = tonic.inhibition,
                                             noise.activation.spread = noise.activation.spread)
  else 
    activation.inhibition.layer <- add.noise(c(n.features, n.locations),
                                             noise.activtion.m = tonic.inhibition,
                                             noise.activation.spread = 0)
  
  weights.inhibition.from.input <- array (generate.random.numbers.around.target (n.features * n.locations * n.features * n.locations,
                                                                                 standard.weight, noise.weight.spread),
                                          dim=c(n.features, n.locations, n.features, n.locations),
                                          dimnames=list(paste("input.feature", 1:n.features, sep="."),
                                                        paste("input.location", 1:n.locations, sep="."),
                                                        paste("inhibition.feature", 1:n.features, sep="."),
                                                        paste("inhibition.location", 1:n.locations, sep="."))) 
  
  # Set connections to same features to a negative value:
  # This is disinhibition
  for (i.feat in 1:n.features){
    weights.inhibition.from.input[i.feat,,i.feat,] <- -weights.inhibition.from.input[i.feat,,i.feat,]    
  }  
  
  # Set self-disinhibition to zero:
  # No self-disinhibition
  for (i.feat in 1:n.features){
    for (i.loc in 1:n.locations){
      weights.inhibition.from.input[i.feat,i.loc,i.feat,i.loc] <- 0
    }
  }
  
  
  weights.inhibition.from.selfdisinhibition <- array (0,
                                                      dim=c(n.features, n.locations, n.features, n.locations),
                                                      dimnames=list(paste("selfdisinhibition.feature", 1:n.features, sep="."),
                                                                    paste("selfdisinhibition.location", 1:n.locations, sep="."),
                                                                    paste("inhibition.feature", 1:n.features, sep="."),
                                                                    paste("inhibition.location", 1:n.locations, sep="."))) 
  # Set connections to same features to a negative value:
  # This is disinhibition
  for (i.feat in 1:n.features){
    weights.inhibition.from.selfdisinhibition[i.feat,,i.feat,] <- -generate.random.numbers.around.target (1, standard.weight, noise.weight.spread)
  }  
  
  # thresholded neurons receiving input from copy layer
  activation.output.layer <- matrix (generate.random.numbers.around.target (n.features * n.locations,
                                                                            0, noise.activation.spread),
                                     nrow = n.features,
                                     ncol = n.locations)
  weights.output.from.copy <- matrix (generate.random.numbers.around.target (n.features * n.locations,
                                                                             standard.weight, noise.weight.spread),
                                      nrow = n.features,
                                      ncol = n.locations)
  
  
  # Now run the simulation
  
  # If the input is just a single matrix (i.e., the simultaneous case, 
  # put it into a list)
  if (!is.list (input.matrix.list)){
    if (is.matrix(input.matrix.list))
      input.matrix.list <- list (input.matrix.list)
    else
      stop ("First argument must be a matrix or a list of matrices")
  }
  
  for (current.input.matrix in input.matrix.list){
    
    # Update self-inhibition here, so it takes the input from the 
    # previous step
    activation.selfdisinhibition.layer <- weights.selfdisinhibition.from.input *  
      activation.input.layer
    if (ADD.NOISE.BEYOND.INPUT)
      activation.selfdisinhibition.layer <- add.noise (activation.selfdisinhibition.layer,
                                                       noise.activation.spread = noise.activation.spread)
    
    # Now present next input
    activation.input.layer <- current.input.matrix
    activation.input.layer <- add.noise (activation.input.layer,
                                         noise.activation.spread = noise.activation.spread)
    
    if (ADD.NOISE.BEYOND.INPUT){
      activation.inhibition.layer <-  add.noise(c(n.features, n.locations),
                                                noise.activtion.m = tonic.inhibition,
                                                noise.activation.spread = noise.activation.spread)
    } else {
      activation.inhibition.layer <-  add.noise(c(n.features, n.locations),
                                                noise.activtion.m = tonic.inhibition,
                                                noise.activation.spread = 0)
    }
    
    # Input to inhibition from input
    for (i.feat in 1:n.features){
      for (i.loc in 1:n.locations){
        activation.inhibition.layer <- activation.inhibition.layer +
          weights.inhibition.from.input[i.feat,i.loc,,] * activation.input.layer[i.feat,i.loc]
      }
    }
    # Input to inhibition from self-disinhibition
    for (i.feat in 1:n.features){
      for (i.loc in 1:n.locations){
        activation.inhibition.layer <- activation.inhibition.layer +
          weights.inhibition.from.selfdisinhibition[i.feat,i.loc,,] * activation.selfdisinhibition.layer[i.feat,i.loc]
      }
    }
    
    activation.inhibition.layer <- restrict.value.range (activation.inhibition.layer)
    
    activation.copy.layer <- (weights.copy.from.input * activation.input.layer) + 
      weights.copy.from.inhibition * activation.inhibition.layer
    
    activation.copy.layer <- restrict.value.range (activation.copy.layer)
    
    activation.output.layer <- 1 * ((weights.output.from.copy * activation.copy.layer) >= threshold.output)
  }
  
  return (list (activation.copy.layer = activation.copy.layer, 
                activation.output.layer = activation.output.layer,
                activation.inhibition.layer = activation.inhibition.layer))
}

#######################
# Run simulations for sequential presentation

# Generate inputs for the noiseless situation
# input.list.repeated <- list (c(1, rep(0, N.FEATURES-1)),
#                              c(1, rep(0, N.FEATURES-1)))
# input.list.repeated <- lapply (input.list.repeated,
#                                function (X) X + generate.random.numbers.around.target (N.FEATURES, 
#                                                        0, NOISE.ACTIVATION.SPREAD))
# 
# input.list.notrepeated <- list (c(1, rep(0, N.FEATURES-1)),
#                                 c(0, 1, rep(0, N.FEATURES-2)))
# input.list.notrepeated <- lapply (input.list.notrepeated,
#                                   function (X) X + generate.random.numbers.around.target (N.FEATURES, 
#                                                           0, NOISE.ACTIVATION.SPREAD))
# run.simulation.sequential (input.list.repeated)
# run.simulation.sequential (input.list.notrepeated) 


# Generate inputs for varying noise
input.list.repeated.no.noise <- list (c(1, rep(0, N.FEATURES-1)),
                                      c(1, rep(0, N.FEATURES-1)))

input.list.notrepeated.no.noise <- list (c(1, rep(0, N.FEATURES-1)),
                                         c(0, 1, rep(0, N.FEATURES-2)))

results.sequential <- data.frame (subj = numeric (), 
                                  noise.spread = numeric(), 
                                  activation.repeated = numeric(), 
                                  activation.not.repeated = numeric(),
                                  accuracy.repeated = numeric(), 
                                  accuracy.not.repeated = numeric())

for (subj in 1:N.SUBJ){
  for (current.noise.spread in NOISE.ACTIVATION.SPREAD.SEQ){
    
    input.list.repeated <- lapply (input.list.repeated.no.noise,
                                   function (X) X + generate.random.numbers.around.target (N.FEATURES, 
                                                                                           0, current.noise.spread))
    
    input.list.notrepeated <- lapply (input.list.notrepeated.no.noise,
                                      function (X) X + generate.random.numbers.around.target (N.FEATURES, 
                                                                                              0, current.noise.spread))
    
    output.sequential.repeated <-run.simulation.sequential (input.list.repeated, 
                                                            noise.activation.spread = current.noise.spread)
    output.sequential.notrepeated <- run.simulation.sequential (input.list.notrepeated,
                                                                noise.activation.spread = current.noise.spread)
    
    results.sequential <- rbind (results.sequential,
                                 cbind (subj = subj,
                                        noise.spread = current.noise.spread,
                                        activation.repeated = extract.target.activation.from.copy.layer(output.sequential.repeated,
                                                                                                        input.list.repeated.no.noise,
                                                                                                        "raw"),
                                        activation.not.repeated = extract.target.activation.from.copy.layer(output.sequential.notrepeated,
                                                                                                            input.list.notrepeated.no.noise,
                                                                                                            "raw"),
                                        accuracy.repeated = extract.target.activation.from.output.layer(output.sequential.repeated,
                                                                                                        input.list.repeated.no.noise),
                                        accuracy.not.repeated = 1 - extract.target.activation.from.output.layer(output.sequential.notrepeated,
                                                                                                                input.list.notrepeated.no.noise)))
  }
}

results.sequential.long <- gather(results.sequential, "sequence_type", activation, 
                                  activation.repeated:accuracy.not.repeated, 
                                  factor_key=TRUE)

plot.activation.sequential <- make.plot (results.sequential.long[-grep("accuracy", results.sequential.long$sequence_type),],
                                         # The rest are default values
                                         measure.var = "activation",
                                         group.vars = c("noise.spread", "sequence_type"),
                                         facet.var = "sequence_type",
                                         x.var = "noise.spread",
                                         xlab = "Noise level", ylab = "Activation (raw)")

plot.accuracy.sequential <- make.plot (results.sequential.long[grep("accuracy", results.sequential.long$sequence_type),],
                                       # The rest are default values
                                       measure.var = "activation",
                                       group.vars = c("noise.spread", "sequence_type"),
                                       facet.var = "sequence_type",
                                       x.var = "noise.spread",
                                       xlab = "Noise level", ylab = "Accuracy")

pdf ("figs/sequential.activation.pdf")
print (plot.activation.sequential)
dev.off ()

pdf ("figs/sequential.accuracy.pdf")
print (plot.accuracy.sequential)
dev.off ()


######################

# Run simulation for simultaneous presentation
# Generate inputs for the noiseless situation
# input.matrix.repeated  <- matrix (generate.random.numbers.around.target (N.FEATURES * N.LOCATIONS, 
#                                          0, NOISE.ACTIVATION.SPREAD),
#                                   nrow = N.FEATURES,
#                                   ncol = N.LOCATIONS)
# input.matrix.repeated[1,1:2] <- 1
# 
# input.matrix.notrepeated  <- matrix (generate.random.numbers.around.target (N.FEATURES * N.LOCATIONS, 
#                                             0, NOISE.ACTIVATION.SPREAD),
#                                      nrow = N.FEATURES,
#                                      ncol = N.LOCATIONS)
# input.matrix.notrepeated[1,1] <- 1
# input.matrix.notrepeated[2,2] <- 1
# 
# run.simulation.simultaneous.topographic.inhibition(input.matrix.repeated)
# run.simulation.simultaneous.topographic.inhibition(input.matrix.notrepeated)   


# Generate inputs for varying noise
results.simultaneous <- data.frame (subj = numeric (), 
                                    noise.spread = numeric(),
                                    activation.repeated = numeric(), 
                                    activation.not.repeated = numeric(),
                                    accuracy.not.repeated = numeric())

for (subj in 1:N.SUBJ){
  for (current.noise.spread in NOISE.ACTIVATION.SPREAD.SEQ){
    
    input.matrix.repeated  <- matrix (generate.random.numbers.around.target (N.FEATURES * N.LOCATIONS, 
                                                                             0, current.noise.spread),
                                      nrow = N.FEATURES,
                                      ncol = N.LOCATIONS)
    input.matrix.repeated[1,1:2] <- 1
    
    input.matrix.notrepeated  <- matrix (generate.random.numbers.around.target (N.FEATURES * N.LOCATIONS, 
                                                                                0, current.noise.spread),
                                         nrow = N.FEATURES,
                                         ncol = N.LOCATIONS)
    input.matrix.notrepeated[1,1] <- 1
    input.matrix.notrepeated[2,2] <- 1
    
    output.simultaneous.repeated <- run.simulation.simultaneous.topographic.inhibition(input.matrix.repeated,
                                                                                       noise.activation.spread = current.noise.spread)
    output.simultaneous.notrepeated <- run.simulation.simultaneous.topographic.inhibition(input.matrix.notrepeated,
                                                                                          noise.activation.spread = current.noise.spread)   
    
    results.simultaneous <- rbind (results.simultaneous,
                                   cbind (subj = subj,
                                          noise.spread = current.noise.spread,
                                          activation.repeated = extract.target.activation.from.copy.layer(output.simultaneous.repeated,
                                                                                                          c(1,2),
                                                                                                          "raw"),
                                          activation.not.repeated = extract.target.activation.from.copy.layer(output.simultaneous.notrepeated,
                                                                                                              c(2,2),
                                                                                                              "raw"),
                                          accuracy.repeated = extract.target.activation.from.output.layer(output.simultaneous.repeated,
                                                                                                          c(1,2)),
                                          accuracy.not.repeated = 1 - extract.target.activation.from.output.layer(output.simultaneous.notrepeated,
                                                                                                                  c(2,2))))
  }
}

results.simultaneous.long <- gather(results.simultaneous, "sequence_type", activation, 
                                    activation.repeated:accuracy.not.repeated, 
                                    factor_key=TRUE)


plot.activation.simultaneous <- make.plot (results.simultaneous.long[-grep("accuracy", results.simultaneous.long$sequence_type),],
                                           # The rest are default values
                                           measure.var = "activation",
                                           group.vars = c("noise.spread", "sequence_type"),
                                           facet.var = "sequence_type",
                                           x.var = "noise.spread",
                                           xlab = "Noise level", ylab = "Activation (raw)")

plot.accuracy.simultaneous <- make.plot (results.simultaneous.long[grep("accuracy", results.simultaneous.long$sequence_type),],
                                         # The rest are default values
                                         measure.var = "activation",
                                         group.vars = c("noise.spread", "sequence_type"),
                                         facet.var = "sequence_type",
                                         x.var = "noise.spread",
                                         xlab = "Noise level", ylab = "Accuracy")


pdf ("figs/simultaneous.activation.pdf")
print (plot.activation.simultaneous)
dev.off ()

pdf ("figs/simultaneous.accuracy.pdf")
print (plot.accuracy.simultaneous)
dev.off ()



#########################

# 3. Run simulations for combined presentation code

# 3.1 Simultaneous case
#run.simulation.combined(input.matrix.repeated)
#run.simulation.combined(input.matrix.notrepeated)   

results.simultaneous.comb <- data.frame (subj = numeric (), 
                                         noise.spread = numeric(), 
                                         activation.repeated.different.location = numeric(), 
                                         activation.not.repeated = numeric(),
                                         accuracy.repeated.different.location = numeric(), 
                                         accuracy.not.repeated = numeric())

for (subj in 1:N.SUBJ){
  for (current.noise.spread in NOISE.ACTIVATION.SPREAD.SEQ){
    
    input.matrix.repeated.comb  <- matrix (generate.random.numbers.around.target (N.FEATURES * N.LOCATIONS, 
                                                                                  0, current.noise.spread),
                                           nrow = N.FEATURES,
                                           ncol = N.LOCATIONS)
    input.matrix.repeated.comb[1,1:2] <- 1
    
    input.matrix.notrepeated.comb  <- matrix (generate.random.numbers.around.target (N.FEATURES * N.LOCATIONS, 
                                                                                     0, current.noise.spread),
                                              nrow = N.FEATURES,
                                              ncol = N.LOCATIONS)
    input.matrix.notrepeated.comb[1,1] <- 1
    input.matrix.notrepeated.comb[2,2] <- 1
    
    output.simultaneous.repeated.comb <- run.simulation.combined(input.matrix.repeated.comb,
                                                                 noise.activation.spread = current.noise.spread)
    output.simultaneous.notrepeated.comb <- run.simulation.combined(input.matrix.notrepeated.comb,
                                                                    noise.activation.spread = current.noise.spread)   
    
    results.simultaneous.comb <- rbind (results.simultaneous.comb,
                                        cbind (subj = subj,
                                               noise.spread = current.noise.spread,
                                               activation.repeated.different.location = extract.target.activation.from.copy.layer(output.simultaneous.repeated.comb,
                                                                                                                                  c(1,2),
                                                                                                                                  "raw"),
                                               activation.not.repeated = extract.target.activation.from.copy.layer(output.simultaneous.notrepeated.comb,
                                                                                                                   c(2,2),
                                                                                                                   "raw"),
                                               accuracy.repeated.different.location = extract.target.activation.from.output.layer(output.simultaneous.repeated.comb,
                                                                                                                                  c(2,2)),
                                               accuracy.not.repeated = 1 - extract.target.activation.from.output.layer(output.simultaneous.notrepeated.comb,
                                                                                                                       c(2,2))))
  }
}

results.simultaneous.comb.long <- gather(results.simultaneous.comb, "sequence_type", activation, 
                                         activation.repeated.different.location:accuracy.not.repeated, 
                                         factor_key=TRUE)


plot.activation.simultaneous.comb <- make.plot (results.simultaneous.comb.long[-grep("accuracy", results.simultaneous.comb.long$sequence_type),],
                                                # The rest are default values
                                                measure.var = "activation",
                                                group.vars = c("noise.spread", "sequence_type"),
                                                facet.var = "sequence_type",
                                                x.var = "noise.spread",
                                                xlab = "Noise level", ylab = "Activation (raw)")

plot.accuracy.simultaneous.comb <- make.plot (results.simultaneous.comb.long[grep("accuracy", results.simultaneous.comb.long$sequence_type),],
                                              # The rest are default values
                                              measure.var = "activation",
                                              group.vars = c("noise.spread", "sequence_type"),
                                              facet.var = "sequence_type",
                                              x.var = "noise.spread",
                                              xlab = "Noise level", ylab = "Accuracy")


pdf ("figs/simultaneous.activation.comb.pdf")
print (plot.activation.simultaneous.comb)
dev.off ()

pdf ("figs/simultaneous.accuracy.comb.pdf")
print (plot.activation.simultaneous.comb)
dev.off ()


####

# 3.2 Sequential case; here we differentiate between repetitions in same and different positions 

input.matrix.base  <- matrix (generate.random.numbers.around.target (N.FEATURES * N.LOCATIONS, 
                                                                     0, NOISE.ACTIVATION.SPREAD),
                              nrow = N.FEATURES,
                              ncol = N.LOCATIONS)

input.matrix.list.notrepeated.no.noise <- list (input.matrix.base, input.matrix.base)
input.matrix.list.notrepeated.no.noise[[1]][1,1] <- 1
input.matrix.list.notrepeated.no.noise[[2]][2,2] <- 1

input.matrix.list.repeated.same.location.no.noise <- list (input.matrix.base, input.matrix.base)
input.matrix.list.repeated.same.location.no.noise[[1]][1,1] <- 1
input.matrix.list.repeated.same.location.no.noise[[2]][1,1] <- 1


input.matrix.list.repeated.different.location.no.noise <- list (input.matrix.base, input.matrix.base)
input.matrix.list.repeated.different.location.no.noise[[1]][1,1] <- 1
input.matrix.list.repeated.different.location.no.noise[[2]][1,2] <- 1

# run.simulation.combined(input.matrix.list.repeated.same.location)
# run.simulation.combined(input.matrix.list.repeated.different.location)
# run.simulation.combined(input.matrix.list.notrepeated)


results.sequential.comb <- data.frame (subj = numeric (), 
                                       noise.spread = numeric(), 
                                       activation.repeated.same.location = numeric(), 
                                       activation.repeated.different.location = numeric(), 
                                       activation.not.repeated = numeric(),
                                       accuracy.repeated.same.location = numeric(),
                                       accuracy.repeated.different.location = numeric(), 
                                       accuracy.not.repeated = numeric())

for (subj in 1:N.SUBJ){
  for (current.noise.spread in NOISE.ACTIVATION.SPREAD.SEQ){
    
    input.list.repeated.same.location.comb <- lapply (input.matrix.list.repeated.same.location.no.noise,
                                                      function (X) X + generate.random.numbers.around.target (N.FEATURES, 
                                                                                                              0, current.noise.spread))
    
    input.list.repeated.different.location.comb <- lapply (input.matrix.list.repeated.different.location.no.noise,
                                                           function (X) X + generate.random.numbers.around.target (N.FEATURES, 
                                                                                                                   0, current.noise.spread))
    
    input.list.notrepeated.comb <- lapply (input.matrix.list.notrepeated.no.noise,
                                           function (X) X + generate.random.numbers.around.target (N.FEATURES, 
                                                                                                   0, current.noise.spread))
    
    output.sequential.repeated.same.location.comb <- run.simulation.combined (input.list.repeated.same.location.comb, 
                                                                              noise.activation.spread = current.noise.spread)
    output.sequential.repeated.different.location.comb <- run.simulation.combined (input.list.repeated.different.location.comb, 
                                                                                   noise.activation.spread = current.noise.spread)
    output.sequential.notrepeated.comb <- run.simulation.combined (input.list.notrepeated.comb,
                                                                   noise.activation.spread = current.noise.spread)
    
    results.sequential.comb <- rbind (results.sequential.comb,
                                      cbind (subj = subj,
                                             noise.spread = current.noise.spread,
                                             activation.repeated.same.location = extract.target.activation.from.copy.layer(output.sequential.repeated.same.location.comb,
                                                                                                                           input.matrix.list.repeated.same.location.no.noise,
                                                                                                                           "raw"),
                                             activation.repeated.different.location = extract.target.activation.from.copy.layer(output.sequential.repeated.different.location.comb,
                                                                                                                                input.matrix.list.repeated.different.location.no.noise,
                                                                                                                                "raw"),
                                             activation.not.repeated = extract.target.activation.from.copy.layer(output.sequential.notrepeated.comb,
                                                                                                                 input.matrix.list.notrepeated.no.noise,
                                                                                                                 "raw"),
                                             accuracy.repeated.same.location = extract.target.activation.from.output.layer(output.sequential.repeated.same.location.comb,
                                                                                                                           input.matrix.list.repeated.same.location.no.noise),
                                             accuracy.repeated.different.location = extract.target.activation.from.output.layer(output.sequential.repeated.different.location.comb,
                                                                                                                                input.matrix.list.repeated.different.location.no.noise),
                                             accuracy.not.repeated = 1 - extract.target.activation.from.output.layer(output.sequential.notrepeated.comb,
                                                                                                                     input.matrix.list.notrepeated.no.noise)))
  }
}

results.sequential.comb.long <- gather(results.sequential.comb, "sequence_type", activation, 
                                       activation.repeated.same.location:accuracy.not.repeated, 
                                       factor_key=TRUE)

plot.activation.sequential.comb <- make.plot (results.sequential.comb.long[-grep("accuracy", results.sequential.comb.long$sequence_type),],
                                              # The rest are default values
                                              measure.var = "activation",
                                              group.vars = c("noise.spread", "sequence_type"),
                                              facet.var = "sequence_type",
                                              x.var = "noise.spread",
                                              xlab = "Noise level", ylab = "Activation (raw)")

plot.accuracy.sequential.comb <- make.plot (results.sequential.comb.long[grep("accuracy", results.sequential.comb.long$sequence_type),],
                                            # The rest are default values
                                            measure.var = "activation",
                                            group.vars = c("noise.spread", "sequence_type"),
                                            facet.var = "sequence_type",
                                            x.var = "noise.spread",
                                            xlab = "Noise level", ylab = "Accuracy")

pdf ("figs/sequential.comb.activation.pdf")
print (plot.activation.sequential.comb)
dev.off ()

pdf ("figs/sequential.comb.accuracy.pdf")
print (plot.accuracy.sequential.comb)
dev.off ()

# These are all plots for convenience :)

# plot.activation.sequential
# plot.activation.sequential.comb ### IS THIS THE CORRECT ONE? 
# plot.activation.simultaneous
# plot.activation.simultaneous.comb
# 
# multiplot (plot.activation.sequential,
#            plot.activation.sequential.comb, 
#            plot.activation.simultaneous,
#            plot.activation.simultaneous.comb,
#            cols = 2)
# 
# 
# plot.accuracy.sequential
# plot.accuracy.sequential.comb 
# plot.accuracy.simultaneous
# plot.accuracy.sequential.comb
# 
# multiplot (plot.accuracy.sequential,
#            plot.accuracy.sequential.comb, 
#            plot.accuracy.simultaneous, 
#            plot.accuracy.sequential.comb,
#            cols =2)

results.overall.long <- rbind (cbind.data.frame(model.type = "dedicated",
                                                sim.type = "sequential",
                                                measure.type = "activation",
                                                results.sequential.long),
                               cbind.data.frame(model.type = "dedicated",
                                                sim.type = "simultaneous",
                                                measure.type = "activation",
                                                results.simultaneous.long),
                               cbind.data.frame(model.type = "combined",
                                                sim.type = "sequential",
                                                measure.type = "activation",
                                                results.sequential.comb.long),
                               cbind.data.frame(model.type = "combined",
                                                sim.type = "simultaneous",
                                                measure.type = "activation",
                                                results.simultaneous.comb.long))

results.overall.long$measure.type <- factor (results.overall.long$measure.type,
                                             levels = c("activation", "accuracy"))
results.overall.long[grep("accuracy",
                          results.overall.long$sequence_type),]$measure.type <- "accuracy"

results.overall.long$measure.type <- mapvalues(results.overall.long$measure.type, 
                                               from = c("activation", "accuracy"),
                                               to = c("activation (copy layer)", "accuracy (output layer)"))

results.overall.long$sequence_type <- factor (results.overall.long$sequence_type,
                                              c("activation.repeated", 
                                                "activation.repeated.same.location",
                                                "activation.repeated.different.location",
                                                "activation.not.repeated",
                                                "accuracy.repeated", 
                                                "accuracy.repeated.same.location",
                                                "accuracy.repeated.different.location",
                                                "accuracy.not.repeated"))


plot.dedicated <- make.plot (results.overall.long[results.overall.long$model.type=="dedicated",],
                             # The rest are default values
                             measure.var = "activation",
                             group.vars = c("noise.spread", "sequence_type"),
                             facet.var = "sequence_type",
                             x.var = "noise.spread",
                             xlab = "Noise level", 
                             facet.grid.vars = c("measure.type", "sim.type"))

pdf ("figs/dedicated.pdf")
print (plot.dedicated)
dev.off ()

plot.combined <- make.plot (results.overall.long[results.overall.long$model.type=="combined",],
                            # The rest are default values
                            measure.var = "activation",
                            group.vars = c("noise.spread", "sequence_type"),
                            facet.var = "sequence_type",
                            x.var = "noise.spread",
                            xlab = "Noise level", 
                            facet.grid.vars = c("measure.type", "sim.type"))

pdf ("figs/combined.pdf")
print (plot.combined)
dev.off ()



plot.activation.only <- make.plot (results.overall.long[results.overall.long$measure.type=="activation (copy layer)",],
                                   # The rest are default values
                                   measure.var = "activation",
                                   group.vars = c("noise.spread", "sequence_type"),
                                   facet.var = "sequence_type",
                                   x.var = "noise.spread",
                                   xlab = "Noise level", 
                                   facet.grid.vars = c("sim.type", "model.type"))

pdf ("figs/activation.4panels.pdf")
print (plot.activation.only)
dev.off ()

plot.activation.only.dedicated <- make.plot (results.overall.long[
  results.overall.long$measure.type=="activation (copy layer)" &
  results.overall.long$model.type=="dedicated",],
                                   # The rest are default values
                                   measure.var = "activation",
                                   group.vars = c("noise.spread", "sequence_type"),
                                   facet.var = "sequence_type",
                                   x.var = "noise.spread",
                                   xlab = "Noise level", 
                                   facet.grid.vars = c("sim.type"))

pdf ("figs/activation.dedicated.pdf")
print (plot.activation.only.dedicated)
dev.off ()


plot.activation.only.combined <- make.plot (results.overall.long[
  results.overall.long$measure.type=="activation (copy layer)" &
    results.overall.long$model.type=="combined",],
  # The rest are default values
  measure.var = "activation",
  group.vars = c("noise.spread", "sequence_type"),
  facet.var = "sequence_type",
  x.var = "noise.spread",
  xlab = "Noise level", 
  facet.grid.vars = c("sim.type"))

pdf ("figs/activation.combined.pdf")
print (plot.activation.only.combined)
dev.off ()
