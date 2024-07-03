# all aux

# Load necessary libraries
library(sp)
library(jpeg)
library(zoo)
library(scales)
library(shiny)

# helper fucntionss


Boundary <- function(X, min, max) {
  # Replace values in X that are less than min or greater than max with NA
  X <- ifelse(X < min | X > max, NA, X)
  return(X)  # Return the modified vector X
}

robust <- function(x, Hz) {
  # Calculate the run-length encoding of the input vector x
  mis <- rle(ifelse(is.na(x), 0, 1))
  
  # Return the mean length of non-missing segments divided by the sampling rate (Hz)
  return(mean(mis$lengths[mis$values == 1]) / Hz)
}

precision <- function(X, Hz) {
  # Calculate the number of windows based on the length of X and the sampling rate
  end <- length(X) / (Hz / 10)
  
  # Create a vector indicating the window number for each data point
  window <- rep(1:end, each = (Hz / 10))
  
  # Apply the MFW function to each window of data points
  Smooth <- as.vector(unlist(by(X[1:length(window)], window, MFW)))
  
  # Calculate and return the mean absolute difference between the original and smoothed data
  return(mean(abs(X[1:length(window)] - Smooth), na.rm = TRUE))
}

Speed_Deg <- function(X, Y, distance, height_mm, width_mm, height_px, width_px, Hz) {
  # Convert horizontal pixel movements to degrees
  hor <- atan((width_mm / 2) / distance) * (180 / pi) * 2 / width_px * X
  
  # Convert vertical pixel movements to degrees
  ver <- atan((height_mm / 2) / distance) * (180 / pi) * 2 / height_px * Y
  
  # Calculate the speed in degrees per second
  speed <- sqrt(diff(hor, 2) ^ 2 + diff(ver, 2) ^ 2) * (Hz / 2)
  
  # Return the speed vector with leading and trailing values to match original length
  return(c(0.001, speed, 0.001))
}


## Interpolation function
Interpolate <- function(X, Y, D, height_mm, width_mm, height_px, width_px, res_x = res_x, res_y = res_y, Hz = Hz, in_thres = in_thres, thres_dur = thres_dur) {
  
  # Calculate speed of eye movements
  s <- Speed(X, Y, D, height_mm, width_mm, height_px, width_px, res_x = res_x, res_y = res_y, Hz)
  
  # Filter out implausibly high speeds
  s <- ifelse(s > 1000, NA, s)
  
  # Check if there are enough local maxima in speed
  if(length(lomax(s)) < 10) {
    # Not enough data for reliable classification, return default values
    return(list('No Return', 'No Return', 'No Return', 'No Return', 'No Return', 'No Return', 'No Return', 'No Return'))
  } else {
    # Calculate velocity threshold using Mould method
    M <- Mould_vel(s, Hz)
    
    # Classify segments as saccades or fixations based on the threshold
    classification <- ifelse(s > M, 'saccade', 'fixation')
    classification[is.na(classification)] <- 'missing'
    
    # Run-length encode the classification
    CL <- rle(classification)
    d <- data.frame(CL$value, CL$length, c(1, cumsum(CL$length)[-length(CL$length)] + 1), cumsum(CL$length))
    names(d) <- c('index', 'dur', 'start', 'end')
    
    dat_x <- X
    dat_y <- Y
    dat_d <- D
    
    # Interpolate missing data points
    for(i in which(d$index == 'missing')) {
      if(i > 1 & i < dim(d)[1] & d[i, 2] < (in_thres * (Hz / 1000))) {
        if(d[i + 1, 1] == 'fixation' & d[i - 1, 1] == 'fixation') {
          ii_s <- d[i - 1, 4]
          ii_e <- d[i + 1, 3]
          speed <- Speed(c(dat_x[ii_s], dat_x[ii_s], dat_x[ii_e]), c(dat_y[ii_s], dat_y[ii_s], dat_y[ii_e]), c(dat_d[ii_s], dat_d[ii_s], dat_d[ii_e]), height_mm, width_mm, height_px, width_px, res_x = res_x, res_y = res_y, Hz)
          if(speed[2] < M) {
            dat_x[d[i, 3] : d[i, 4]] <- dat_x[ii_s]
            dat_y[d[i, 3] : d[i, 4]] <- dat_y[ii_s]
            dat_d[d[i, 3] : d[i, 4]] <- dat_d[ii_s]
          }
        }
      }
    }
    
    # Recalculate speed after interpolation
    s <- Speed(dat_x, dat_y, dat_d, height_mm, width_mm, height_px, width_px, res_x = res_x, res_y = res_y, Hz)
    s <- ifelse(s > 1000, NA, s)
    
    # Reclassify segments based on the updated speed
    classification <- ifelse(s > M, 'saccade', 'fixation')
    classification[is.na(classification)] <- 'missing'
    CL <- rle(classification)
    
    # Calculate point of gaze (POG)
    index <- rep.int(1:length(CL$value), CL$lengths)
    POG <- sapply(unique(index[!is.na(index)]), function(i) mean(dist(cbind(dat_x[index == i], dat_y[index == i])), na.rm = TRUE))
    POG[is.na(POG)] <- 0
    mean_x <- as.vector(by(dat_x, index, function(i) mean(i, na.rm = TRUE)))
    mean_y <- as.vector(by(dat_y, index, function(i) mean(i, na.rm = TRUE)))
    
    # Create data frame for fixations and saccades
    d <- data.frame(CL$value, CL$length, c(1, cumsum(CL$length)[-length(CL$length)] + 1), cumsum(CL$length), POG, mean_x, mean_y)
    names(d)[1:4] <- c('index', 'dur', 'start', 'end')
    
    ## Combine fixations 
    dimd_new <- dim(d)[1] + 1
    while(dimd_new != dim(d)[1]){
      dimd_new <- dim(d)[1]
      ## Combine fixations
      classif <- comhull(d, classification, dat_x, dat_y, in_thres, Hz, M = sqrt(M)/10, mean(dat_d, na.rm = TRUE), res_x = res_x, width_mm = width_mm)
      
      CL <- rle(classif[[1]])
      classification <- classif[[1]]
      index <- rep.int(1:length(CL$value), CL$lengths)
      dat_x <- classif[[2]]
      dat_y <- classif[[3]]
      POG <- sapply(unique(index[!is.na(index)]), function(i) mean(dist(cbind(dat_x[index == i], dat_y[index == i])), na.rm = TRUE))
      POG[is.na(POG)] <- 0
      mean_x <- as.vector(by(dat_x, index, function(i) mean(i, na.rm = TRUE)))
      mean_y <- as.vector(by(dat_y, index, function(i) mean(i, na.rm = TRUE)))
      
      d <- data.frame(CL$value, CL$length, c(1, cumsum(CL$length)[-length(CL$length)] + 1), cumsum(CL$length), POG, mean_x, mean_y)
      names(d)[1:4] <- c('index', 'dur', 'start', 'end')
    }
    
    # Remove short fixations and check for combinations again
    for(i in which(CL$value == 'fixation' & CL$length < (Hz / 1000 * thres_dur))){
      classification[((cumsum(CL$length) - CL$length) + 1)[i] : cumsum(CL$length)[i]] <- 'saccade' 
    }
    
    CL <- rle(classification)
    index <- rep.int(1:length(CL$value), CL$lengths)
    POG <- sapply(unique(index[!is.na(index)]), function(i) mean(dist(cbind(dat_x[index == i], dat_y[index == i])), na.rm = TRUE))
    POG[is.na(POG)] <- 0
    mean_x <- as.vector(by(dat_x, index, function(i) mean(i, na.rm = TRUE)))
    mean_y <- as.vector(by(dat_y, index, function(i) mean(i, na.rm = TRUE)))
    
    d <- data.frame(CL$value, CL$length, c(1, cumsum(CL$length)[-length(CL$length)] + 1), cumsum(CL$length), POG, mean_x, mean_y)
    names(d)[1:4] <- c('index', 'dur', 'start', 'end')
    
    classification <- comhull(d, classification, dat_x, dat_y, in_thres, Hz, M = sqrt(M)/10, D = mean(dat_d, na.rm = TRUE), res_x = res_x, width_mm = width_mm)
    
    clas <- classification[[1]]
    CL <- rle(clas)
    dat_x <- classification[[2]]
    dat_y <- classification[[3]]
    
    index <- rep.int(1:length(CL$value), CL$lengths)
    mean_x <- as.vector(by(dat_x, index, function(i) mean(i, na.rm = TRUE)))
    mean_y <- as.vector(by(dat_y, index, function(i) mean(i, na.rm = TRUE)))
    
    index <- CL$value
    end <- cumsum(CL$length) * (1000 / Hz)
    dur <- CL$length * (1000 / Hz)
    start <- (end - dur) + 1
    
    d <- data.frame(index, dur, start, end, mean_x, mean_y)
    d <- data.frame(d, order = 1:dim(d)[1])
    
    return(list(dat_x, dat_y, dat_d, d, M, s, clas, 'Return'))
  }
}



posthocCheck <- function(classification, x, y) {
  # Run-length encode the classification vector
  class <- rle(classification)
  
  # Create a data frame summarizing the run-length encoding
  simple <- data.frame(class$values, class$lengths, 
                       c(1, cumsum(class$lengths) + 1)[-(length(class$values) + 1)], 
                       cumsum(class$lengths))
  
  # Number of fixations
  lf <- length(which(class$values == 'f'))
  
  # If there are any fixations
  if(lf > 0) {
    # Initialize vectors for start and end positions and mean positions
    x_start <- y_start <- x_end <- y_end <- mean_x <- mean_y <- POGvar_px <- numeric()
    
    # Iterate through each segment to calculate start, end, and mean positions
    for(i in 1:dim(simple)[1]) {
      x_start <- c(x_start, x[simple[i, 3]])
      y_start <- c(y_start, y[simple[i, 3]])
      x_end <- c(x_end, x[simple[i, 4]])
      y_end <- c(y_end, y[simple[i, 4]])
      mean_x <- c(mean_x, mean(x[simple[i, 3] : simple[i, 4]]))
      mean_y <- c(mean_y, mean(y[simple[i, 3] : simple[i, 4]]))
      POGvar_px <- c(POGvar_px, mean(dist(cbind(x[simple[i, 3] : simple[i, 4]], y[simple[i, 3] : simple[i, 4]]))))
    }
    
    # Combine the calculated values into the summary data frame
    simple <- cbind(simple, x_start, y_start, x_end, y_end, mean_x, mean_y, POGvar_px)
  }
  
  ## Check for misclassifications of fixations
  if(length(which(simple[,1] == 'f')) > 1) {
    # Calculate distances between consecutive fixations
    dis <- dist(simple[which(simple[,1] == 'f'), 9:10])
    
    # Extract diagonal elements to get distances between consecutive fixations
    if(length(dis) == 1) {
      dis.between <- dis
    } else {
      dis.between <- diag(as.matrix(dis)[-1, -dim(as.matrix(dis))])
    }
    
    # Get within-fixation point of gaze variability
    dis.within <- simple[which(simple[,1] == 'f'), 11]
    
    # Iterate through distances between fixations
    for(i in 1:length(dis.between)) {
      # If the distance between fixations is less than the combined within-fixation variability
      if(dis.between[i] < dis.within[i] + dis.within[i + 1]) {
        index <- c(simple[which(simple[,1] == 'f')[i], 3], simple[which(simple[,1] == 'f')[i + 1], 4])
        classification[index[1] : index[2]] <- 'f'
        
        # Replace missing x and y values with the mean of the segment
        x[index[1] : index[2]][is.na(x[index[1] : index[2]])] <- mean(x[index[1] : index[2]], na.rm = TRUE)
        y[index[1] : index[2]][is.na(y[index[1] : index[2]])] <- mean(y[index[1] : index[2]], na.rm = TRUE)
      }
    }
  }
  
  return(list(classification, x, y))
}


simplify <- function(classification, x, y, Hz, D, width_px, width_mm, extra, extra_var) {
  # Run-length encode the classification vector
  class <- rle(classification)
  
  # Create a data frame summarizing the run-length encoding
  simple <- data.frame(class$values, class$lengths, 
                       c(1, cumsum(class$lengths) + 1)[-(length(class$values) + 1)], 
                       cumsum(class$lengths))
  
  # If there are any fixations
  if(length(which(class$values == 'f')) > 0) {
    # Initialize vectors for start and end positions, mean positions, POG variability, and RMS
    x_start <- y_start <- x_end <- y_end <- mean_x <- mean_y <- POGvar <- RMS <- numeric()
    
    # Iterate through each segment to calculate start, end, and mean positions
    for(i in 1:dim(simple)[1]) {
      x_start <- c(x_start, x[simple[i, 3]])
      y_start <- c(y_start, y[simple[i, 3]])
      x_end <- c(x_end, x[simple[i, 4]])
      y_end <- c(y_end, y[simple[i, 4]])
      mean_x <- c(mean_x, mean(x[simple[i, 3] : simple[i, 4]]))
      mean_y <- c(mean_y, mean(y[simple[i, 3] : simple[i, 4]]))
      m <- as.matrix(dist(cbind(c(mean_x[length(mean_x)], x[simple[i, 3] : simple[i, 4]]), c(mean_y[length(mean_y)], y[simple[i, 3] : simple[i, 4]]))))
      RMS <- c(RMS, sqrt(mean((atan((diag(m[-1, -c(1, 2)]) / 2) / mean(D, na.rm = TRUE)) * (180 / pi) * (width_mm / width_px) * 2) ** 2)))
      POGvar <- c(POGvar, mean(m[-1, 1]))
    }
    
    ## Calculate saccade amplitude and transform POGvar from pixels to degrees of visual angle and to standard deviation
    ss <- which(class$values == 's')
    POGvar[ss] <- sqrt((x_start[ss] - x_end[ss]) ^ 2 + (y_start[ss] - y_end[ss]) ^ 2)
    POGsdSacAmp <- atan((POGvar / 2) / mean(D, na.rm = TRUE)) * (180 / pi) * (width_mm / width_px) * 2
    POGsdSacAmp[!ss] <- sqrt(POGsdSacAmp)
    RMS[ss] <- NA
    
    # Update the simple data frame with calculated metrics
    simple <- data.frame(class$values, class$lengths * (1000 / Hz), 
                         c(1, cumsum(class$lengths * (1000 / Hz)) + 1)[-(length(class$values) + 1)], 
                         cumsum(class$lengths * (1000 / Hz)), x_start, y_start, x_end, y_end, mean_x, mean_y, POGsdSacAmp, RMS)
    names(simple)[1:4] <- c('Value', 'Dur', 'Start', 'End')
    
    # Include additional variables if specified
    if(!is.null(extra_var)) {
      for(i in 1:length(extra_var)) {
        simple <- data.frame(simple, extra[i])
        names(simple)[dim(simple)[2]] <- extra_var[i]
      }
    }
  }
  
  # Remove rows with NA values in the first column
  if(length(which(is.na(simple[, 1]))) != 0) {
    simple <- simple[-which(is.na(simple[, 1])), ]
  }
  
  return(simple)
}


#

pnt.in.poly <- function(points, poly) {
  # Ensure the input points and polygon are in the right format
  points <- as.matrix(points)
  poly <- as.matrix(poly)
  
  # Create the polygon object
  polygon <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(poly)), ID = 1)))
  
  # Create the points object
  points <- sp::SpatialPoints(points)
  
  # Check which points are inside the polygon
  inside <- sp::over(points, polygon, returnList = TRUE)
  
  # Convert the result to a logical vector
  inside <- sapply(inside, function(x) !is.null(x))
  
  return(inside)
}








