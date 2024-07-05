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


# ORIOGINAL;

simplify <-
  function(classification, x, y, Hz, D, width_px, width_mm, extra, extra_var){
    class <- rle(classification)
    simple <- data.frame(class$values, class$lengths, 
                         c(1, cumsum(class$lengths) + 1)[-(length(class$values) + 1)], 
                         cumsum(class$lengths))
    if(length(which(class$values == 'f')) > 0){
      x_start <- y_start <- x_end <- y_end <- mean_x <- mean_y <- POGvar <- RMS <- numeric()
      for(i in 1:dim(simple)[1]){
        x_start <- c(x_start, x[simple[i, 3]])
        y_start <- c(y_start, y[simple[i, 3]])
        x_end <- c(x_end, x[simple[i, 4]])
        y_end <- c(y_end, y[simple[i, 4]])
        mean_x <- c(mean_x, mean(x[simple[i,3] : simple[i,4]]))
        mean_y <- c(mean_y, mean(y[simple[i,3] : simple[i,4]]))
        m <- as.matrix(dist(cbind(c(mean_x[length(mean_x)], x[simple[i,3] : simple[i,4]]), c(mean_y[length(mean_y)], y[simple[i,3] : simple[i,4]]))))
        RMS <- c(RMS, sqrt(mean((atan((diag(m[-1,-c(1,2)]) / 2) / mean(D, na.rm = T)) * (180 / pi) * (width_mm / width_px) * 2)**2)))
        POGvar <- c(POGvar, mean(m[-1,1]))
      }
      ## Calculate saccade amplitude and transform POGvar from pixels to degrees of visual angle and to sd
      ss <- which(class$values == 's')
      POGvar[ss] <- sqrt((x_start[ss] - x_end[ss]) ^ 2 + (y_start[ss] - y_end[ss]) ^ 2)
      POGsdSacAmp <- atan((POGvar / 2) / mean(D, na.rm = T)) * (180 / pi) * (width_mm / width_px) * 2
      POGsdSacAmp[!ss] <- sqrt(POGsdSacAmp)
      RMS[ss] <- NA
      
      simple <- data.frame(class$values, class$lengths * (1000/Hz), 
                           c(1, cumsum(class$lengths * (1000/Hz)) + 1)[-(length(class$values) + 1)], 
                           cumsum(class$lengths * (1000/Hz)), x_start, y_start, x_end, y_end, mean_x, mean_y, POGsdSacAmp, RMS)
      names(simple)[1:4] <- c('Value', 'Dur', 'Start', 'End')
      if(!is.null(extra_var)){
        for(i in 1:length(extra_var)){
          simple <- data.frame(simple, extra[i])
          names(simple)[dim(simple)[2]] <- extra_var[i]
        }
      }
    }
    # Remove NA values
    if(length(which(is.na(simple[,1]))) != 0){
      simple <- simple[-which(is.na(simple[,1])),]
    }
    return(simple)
  }


comhull<-
function(d, classification, dat_x, dat_y, in_thres, Hz = Hz, M, D, res_x = res_x, width_mm = width_mm){
  d <- d[d$dur > 1,]
  fix <- tail(which(d$index == 'fixation'), 1)
  count <- length(which(d$index == 'fixation')) - 1
  while(count >= 1){
    fix2 <- which(d$index == 'fixation')[count]
    cvhull <- chull(cbind(dat_x, dat_y)[d[fix,3] : d[fix,4],])
    POLY_FIX <- cbind(dat_x[d[fix,3] : d[fix,4]][cvhull], dat_y[d[fix,3] : d[fix,4]][cvhull])
    
    # Debug print statements
    # print("Dimensions of points:")
    # print(dim(cbind(dat_x, dat_y)[d[fix2,3] : d[fix2,4],]))
    # print("Dimensions of polygon:")
    # print(dim(POLY_FIX))
    
    # Modified line to use the updated pnt.in.poly function
    PNT <- sum(pnt.in.poly(cbind(dat_x, dat_y)[d[fix2,3] : d[fix2,4],], POLY_FIX), na.rm = TRUE)
    
    # print("PNT value:")
    # print(PNT)
    # 
    ## If fixations have the same location, are within interpolation limit and below distance limit combine them.
    if(!is.na(PNT) && PNT != 0){
      dis <- dist(rbind(t(apply(cbind(dat_x, dat_y)[d[fix,3] : d[fix,4],], 2, mean)),
                        t(apply(cbind(dat_x, dat_y)[d[fix2,3] : d[fix2,4],], 2, mean))))
      thres_d <- atan((width_mm/2)/D) * (180/pi) * 2 *(dis/res_x)
      if(thres_d < M & (d[fix,3] - d[fix2,4]) < in_thres * (Hz / 1000)){
        classification[d[fix2,3] : d[fix,4]] <- 'fixation'
        CL <- rle(classification)
        index <- rep.int(1:length(CL$value), CL$lengths)
        POG <- sapply(unique(index[!is.na(index)]), function(i) mean(dist(cbind(dat_x[index == i], dat_y[index == i])), na.rm = T))
        POG[is.na(POG)] <- 0
        mean_x <- as.vector(by(dat_x, index, function(i) mean(i, na.rm = T)))
        mean_y <- as.vector(by(dat_y, index, function(i) mean(i, na.rm = T)))
        
        dat_x[d[fix2,3] : d[fix,4]] <- na.approx(dat_x[d[fix2,3] : d[fix,4]])
        dat_y[d[fix2,3] : d[fix,4]] <- na.approx(dat_y[d[fix2,3] : d[fix,4]])
        
        d <- data.frame(CL$value, CL$length, c(1, cumsum(CL$length)[-length(CL$length)] + 1), cumsum(CL$length), POG, mean_x, mean_y)
        names(d)[1:4] <- c('index', 'dur', 'start', 'end')
        d <- d[d$dur > 1,]
      } 
    }
    fix <- fix2
    count <- count - 1
  }
  return(list(classification, dat_x, dat_y))
}

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


## This function calculates the smooth for the precision function
MFW <- function(x){
  if(sum(is.na(x)) / length(x) > .5){
    return(rep(NA, length(x)))
  } else {
    return(rep(median(x, na.rm = T), length(x)))
  }
}

Speed <-
  function(X, Y, distance, height_mm, width_mm, height_px, width_px, res_x = 1280, res_y = 1024, Hz){
    d2p_px <- sqrt(abs(X - (res_x / 2))**2 + abs(Y - (res_y / 2))**2 + (distance * width_px / width_mm)**2)
    dbp_px <- c(.001, sqrt(diff(X, 2)**2 + diff(Y, 2)**2), .001)
    speed <- atan((dbp_px / 2) / d2p_px) * (180 / pi) * 2 * (Hz / 2)
  }

lomax <-
  function(x){
    which(diff(c(TRUE, diff(x) >= 0, FALSE)) < 0)
  }



Mould_vel <-
  function(speed, Hz, plot = F){
    lmax <- speed[lomax(speed)]
    if(length(lmax) < 10){
      return(NA)
      warning('There are no enough data points to estimate a velocity threshold')
    } else {
      thresholds <- seq(min(lmax), max(lmax), length.out = Hz)
      set <- sapply(thresholds, function(x) {length(which(lmax > x))})
      uni <- seq(length(lmax), 0, length.out = Hz)
      gap <- uni - set
      
      if(Hz < 250) {h <- .1} else {h <- .05}
      
      gap <- predict(loess( gap ~ log(thresholds) ,span = h))
      while (length(lomax(gap)) > 1 & h < 1) {
        h <- h + .01
        gap <- predict(loess( (uni - set) ~ log(thresholds) ,span = h,  surface = 'direct', cell = 1))
      }
      
      if(plot == T) plotMould(uni, set, gap, thresholds, lmax, Hz)
      if(h != 1) return(thresholds[which.max(gap)]) else return(NA)
    }
  }



######summary fuxcntion

# summary.gazepath <-
#   function(object, ..., complete_only = FALSE, fixations_only = FALSE){
#     output <- numeric()
#     end <- dim(object[[16]][[1]])[2]
#     for(i in 1:length(object[[16]])){
#       if(end == 4) end <- dim(object[[16]][[i]])[2]
#       sim <- object[[16]][[i]]
#       l <- length(which(sim[,1] == 'f'))
#       if(l != 0){
#         if(complete_only == TRUE){
#           if(fixations_only == TRUE){
#             index <- complete(sim, 'f')
#             if(length(index) != 0){
#               output <- rbind(output, cbind(sim[index, c(1:4, 9:end)], 1:length(index), i))
#             }
#           } else {
#             if(length(which(sim[,1] == 's')) != 0){
#               index <- sort(c(complete(sim, 'f'), complete(sim, 's')))
#               if(length(index) != 0){
#                 output <- rbind(output, cbind(sim[index, c(1:4, 9:end)], 1:length(index), i))
#               }
#             } 
#           }
#         } else {
#           if(fixations_only == TRUE){
#             output <- rbind(output, cbind(sim[sim[,1] == 'f', c(1:4, 9:end)], 1:l, i))
#           } else {
#             l <- sum(sim[,1] == 'f' | sim[,1] == 's')
#             output <- rbind(output, cbind(sim[sim[,1] == 'f' | sim[,1] == 's',c(1:4, 9:end)], 1:l, i))
#           }
#         }
#       }
#     }
#     if(length(output) == 0){
#       print('There were no fixations or saccades classified, probably data quality of this particpant is very low')
#     } else {
#       names(output)[c(1:4,(end - 3):(end - 2))] <- c('Value', 'Duration', 'Start', 'End', 'Order', 'Trial')
#       row.names(output) <- 1:dim(output)[1]
#       return(output)
#     }
#   }

# SYUMAMRYU THAT INDICATES TRUIAKLS WITH NO FIXATIONS
summary.gazepath <- function(object, ..., complete_only = FALSE, fixations_only = FALSE){
  output <- numeric()
  end <- dim(object[[16]][[1]])[2]
  no_fixation_trials <- c()
  
  for(i in 1:length(object[[16]])){
    if(end == 4) end <- dim(object[[16]][[i]])[2]
    sim <- object[[16]][[i]]
    l <- length(which(sim[,1] == 'f'))
    
    if(l == 0){
      no_fixation_trials <- c(no_fixation_trials, i)
    } else {
      if(complete_only == TRUE){
        if(fixations_only == TRUE){
          index <- complete(sim, 'f')
          if(length(index) != 0){
            output <- rbind(output, cbind(sim[index, c(1:4, 9:end)], 1:length(index), i))
          }
        } else {
          if(length(which(sim[,1] == 's')) != 0){
            index <- sort(c(complete(sim, 'f'), complete(sim, 's')))
            if(length(index) != 0){
              output <- rbind(output, cbind(sim[index, c(1:4, 9:end)], 1:length(index), i))
            }
          } 
        }
      } else {
        if(fixations_only == TRUE){
          output <- rbind(output, cbind(sim[sim[,1] == 'f', c(1:4, 9:end)], 1:l, i))
        } else {
          l <- sum(sim[,1] == 'f' | sim[,1] == 's')
          output <- rbind(output, cbind(sim[sim[,1] == 'f' | sim[,1] == 's',c(1:4, 9:end)], 1:l, i))
        }
      }
    }
  }
  
  if(length(output) == 0){
    print('There were no fixations or saccades classified, probably data quality of this participant is very low')
  } else {
    names(output)[c(1:4,(end - 3):(end - 2))] <- c('Value', 'Duration', 'Start', 'End', 'Order', 'Trial')
    row.names(output) <- 1:dim(output)[1]
    
    # Print trials with no fixations
    if(length(no_fixation_trials) > 0) {
      cat("Trials with no detectable fixations:", paste(no_fixation_trials, collapse = ", "), "\n")
    } else {
      cat("All trials contained detectable fixations.\n")
    }
    
    return(output)
  }
}


# plotting

plot.gazepath <-
  function(x, ..., trial_index = 1){
    i <- trial_index
    object <- x
    if(dim(object[[16]][[i]])[2] == 4){
      warning('There is not enough data to identify fixations and saccades in this trial')
    } else {
      layout(matrix(c(1, 1:3), 2, 2))
      plot(object[[2]][[i]], object[[3]][[i]], xlab = "X", ylab = 'Y', type = 'l', las = 1, ylim = c(max(object[[3]][[i]], na.rm = T), min(object[[3]][[i]], na.rm = T)))
      sim <- object[[16]][[i]]
      points(sim[sim[,1] == 'f', 9:10], pch = letters, cex = 3, col = 4)
      
      plot(object[[2]][[i]], ylim = c(-50, max(object[[14]][[i]], object[[12]][[i]]) + 100), type = 'l', xlab = 'Time (msec)', ylab = 'position', las = 1, xaxt = 'n', col = 5)
      lines(object[[3]][[i]], col = 4)
      axis(1, at = seq(0, length(object[[2]][[i]]), length.out = 6), labels = round(seq(0, length(object[[2]][[i]]) * (1000 / object[[10]]), length.out = 6)))
      fix <- sim[sim[,1] == 'f',3:4] / (1000 / object[[10]])
      rect(fix[,1], -50, fix[,2], 0, col = 3)
      legend('topleft', c('X-coordinates', 'Y-coordinates', 'Fixations'), col = 5:3, lwd = 2, bty = 'n', horiz = TRUE)
      
      if(object[[4]] != 'Tobii' & object[[4]] != 'Eyelink'){
        plot(object[[9]][[i]], typ = 'l', xlab = 'Time (msec)', ylab = 'Speed (deg/s)', las = 1, xaxt = 'n')
        axis(1, at = seq(0, length(object[[9]][[i]]), length.out = 6), labels = round(seq(0, length(object[[9]][[i]]) * (1000 / object[[10]]), length.out = 6)))
        if(object[[4]] == 'Mould.all' | object[[4]] == 'Mould.allDur'){
          segments(0, object[[7]], length(object[[9]][[i]]), object[[7]], col = 2, lwd= 2)
        } else {
          segments(0, object[[7]][[i]], length(object[[9]][[i]]), object[[7]][[i]], col = 2, lwd= 2)
        }
      }
    }
    layout(1)
  }




## Combine succesive fixations based on region, overlapping fixations are combined


#######

comhull <- function(d, classification, dat_x, dat_y, in_thres, Hz = Hz, M, D, res_x = res_x, width_mm = width_mm){
  d <- d[d$dur > 1,]
  fix <- tail(which(d$index == 'fixation'), 1)
  count <- length(which(d$index == 'fixation')) - 1
  while(count >= 1){
    fix2 <- which(d$index == 'fixation')[count]
    cvhull <- chull(cbind(dat_x, dat_y)[d[fix,3] : d[fix,4],])
    POLY_FIX <- cbind(dat_x[d[fix,3] : d[fix,4]][cvhull], dat_y[d[fix,3] : d[fix,4]][cvhull])
    
    # Debug print statements
    # print("Dimensions of points:")
    # print(dim(cbind(dat_x, dat_y)[d[fix2,3] : d[fix2,4],]))
    # print("Dimensions of polygon:")
    # print(dim(POLY_FIX))
    
    # Modified line to use the updated pnt.in.poly function
    PNT <- sum(pnt.in.poly(cbind(dat_x, dat_y)[d[fix2,3] : d[fix2,4],], POLY_FIX), na.rm = TRUE)
    
    # print("PNT value:")
    # print(PNT)
    # 
    ## If fixations have the same location, are within interpolation limit and below distance limit combine them.
    if(!is.na(PNT) && PNT != 0){
      dis <- dist(rbind(t(apply(cbind(dat_x, dat_y)[d[fix,3] : d[fix,4],], 2, mean)),
                        t(apply(cbind(dat_x, dat_y)[d[fix2,3] : d[fix2,4],], 2, mean))))
      thres_d <- atan((width_mm/2)/D) * (180/pi) * 2 *(dis/res_x)
      if(thres_d < M & (d[fix,3] - d[fix2,4]) < in_thres * (Hz / 1000)){
        classification[d[fix2,3] : d[fix,4]] <- 'fixation'
        CL <- rle(classification)
        index <- rep.int(1:length(CL$value), CL$lengths)
        POG <- sapply(unique(index[!is.na(index)]), function(i) mean(dist(cbind(dat_x[index == i], dat_y[index == i])), na.rm = T))
        POG[is.na(POG)] <- 0
        mean_x <- as.vector(by(dat_x, index, function(i) mean(i, na.rm = T)))
        mean_y <- as.vector(by(dat_y, index, function(i) mean(i, na.rm = T)))
        
        dat_x[d[fix2,3] : d[fix,4]] <- na.approx(dat_x[d[fix2,3] : d[fix,4]])
        dat_y[d[fix2,3] : d[fix,4]] <- na.approx(dat_y[d[fix2,3] : d[fix,4]])
        
        d <- data.frame(CL$value, CL$length, c(1, cumsum(CL$length)[-length(CL$length)] + 1), cumsum(CL$length), POG, mean_x, mean_y)
        names(d)[1:4] <- c('index', 'dur', 'start', 'end')
        d <- d[d$dur > 1,]
      } 
    }
    fix <- fix2
    count <- count - 1
  }
  return(list(classification, dat_x, dat_y))
}


pnt.in.poly <- function(points, poly) {
  # Ensure the input points and polygon are in the right format
  points <- as.matrix(points)
  poly <- as.matrix(poly)
  
  # Create the polygon object
  polygon <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(poly)), ID = 1)))
  
  # Create the points object
  points <- sp::SpatialPoints(points)
  
  # Check which points are inside the polygon
  inside <- sp::over(points, polygon, returnList = FALSE)
  
  # Convert the result to a matrix with one column
  inside <- matrix(as.numeric(inside), ncol = 1)
  
  return(inside)
}



