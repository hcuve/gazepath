# Gazepath updated functions - July 8th 2024
This repository is a fork that recreates and updates the gazepath functions for parsing gaze data into fixations and saccades, as the original gazepath package is no longer available and does not run on newer versions of R. I have adapted and updated the functions to run directly in R without the need for the package.


# Example of Usage
run the code u the scripts ```helperfunctions.R``` and  ```GazePathHC.R``` to recreate the necessary fucntions


```
library(tidyverse)

# load the data
#transform the data into  matrix if needed
dta_4gazepath_stim <- as.data.frame(as.matrix(dta_4gazepath_stim))
#transform specific columns that the gazpearth will need into the correct formats
dta_4gazepath_stim[, c(17, 18, 19, 20, 23, 28)] <- lapply(dta_4gazepath_stim[, c(17, 18, 19, 20, 23, 28)], as.numeric)

# Split the data by participant
tmp_stim_by_pp <- split(dta_4gazepath_stim, dta_4gazepath_stim$participant)

# Initialize an empty list to store gazepath results
gazepath_results <- list()

# Initialize an empty dataframe to store summary results
gp_summary_results <- data.frame()

# Loop through each participant's data
for (participant_id in names(tmp_stim_by_pp)) {
  # Get the data for the current participant
  dta <- tmp_stim_by_pp[[participant_id]]

  # Check if the trial column exists and has values
  if (!"trial_no" %in% colnames(dta) || all(is.na(dta$trial_no))) {
    warning(paste("Skipping participant", participant_id, "due to missing or empty 'trial' column"))
    next
  }

  # Ensure the trial column is properly formatted
  dta$trial_no <- as.factor(as.numeric(as.character(dta$trial_no)))

  # Run the gazepath analysis (you can run just just thes ection below if you are running only a single file/particiapnt )
  #
  
  rslt_gp_stim <- gazepath(dta, 
                           x1 = 17, y1 = 19, 
                           x2 = 18, y2 = 20,
                           d1 = 23, d2 = 23,
                           trial = 28, # Ensure this is a character or factor
                           height_px = 1080, height_mm = 295,
                           width_px = 1920, width_mm = 525,
                           extra_var = c("text", "trial_no", "trial_no_gpid_stim", "stimnopath_no_na", "participant"),
                           method = 'gazepath',
                           samplerate = 150)

  # Append the gazepath result to the list
  gazepath_results[[participant_id]] <- rslt_gp_stim

  # Run the summary analysis
  rslt_fix_sacc_gp <- summary.gazepath(rslt_gp_stim)

  # Bind the summary results into the dataframe
  gp_summary_results <- bind_rows(gp_summary_results, rslt_fix_sacc_gp)
}


# Save gazepath results to an RDS file
saveRDS(gazepath_results, "gazepath_results.rds")
saveRDS(gp_summary_results, "gp_summary_results.rds")
saveRDS(tmp_stim_by_pp, "tmp_stim_by_pp.rds")

# Example plots that compare before and after to make sure that parsed data matched the raw data
tmp_stim_by_pp$`140` %>%
  subset(trial_no == 6) %>%
  ggplot(aes(as.numeric(gaze_x_cor_pix), as.numeric(gaze_y_cor_pix))) +
  geom_path() +
  geom_path(aes(y = as.numeric(gaze_y_cor_pix)))+
  scale_y_reverse

tmp_stim_by_pp$`140` %>%
  subset(trial_no == 6) %>%
  ggplot(aes(as.numeric(timerezero), as.numeric(gaze_x_cor_pix))) +
  geom_line() +
  geom_line(aes(y = as.numeric(gaze_y_cor_pix)))

plot.gazepath(gazepath_results$`140`, trial_index = 6)

```


 Arguments
data
The dataframe with at least the raw x- and y-coordinates, the distance to the screen in mm and a trial index.

x1
The column name (between quotes, e.g. 'x1') or the number of the column in the dataframe containing the x-coordinates

y1
The column name (between quotes, e.g. 'y1') or the number of the column in the dataframe containing the y-coordinates

x2
When tracking was binocular, the column name (between quotes, e.g. 'x2') or number of the dataframe containing the x-coordinates of the second eye

y2
When tracking was binocular, the column name (between quotes, e.g. 'y2') or number of the dataframe containing the y-coordinates of the second eye

d1
The column name (between quotes, e.g. 'd2') or numberof the dataframe containing the distance in mm

d2
When tracking was binocular, the column name (between quotes, e.g. 'd2') or number of the dataframe containing the distance in mm of the second eye

trial
The column name (between quotes, e.g. 'TRIAL_INDEX') or number of the dataframe containing the trial or stimuli index

height_px
The height of the stimuli in pixels, can be a single value or a vector of length number of trials when stimuli differ in size per trial

height_mm
The height of the stimuli in mm, can be a single value or a vector of length number of trials when stimuli differ in size per trial

width_px
The width of the stimuli in pixels, can be a single value or a vector of length number of trials when stimuli differ in size per trial

width_mm
The height of the stimuli in pixels, can be a single value or a vector of length number of trials when stimuli differ in size per trial trials

# extra_var A vector of names of the variables that must return in the output file, for example, condition, stimuli name, etc.

The vertical resolution of the monitor in pixels

samplerate
method


#errors
if you are getting errors (e.g. rel related) , it likely has somethign to do with formatting issues of the data, so check datframe, column types and if needed re-do them

# To Do
- Organise the code and fix typos
- Split helper functions into individual scripts
- Recreate the Shiny app


# citation

If you use this, please check and cite the original authors:
van Renswoude, D. R., Raijmakers, M. E., Koornneef, A., Johnson, S. P., Hunnius, S., & Visser, I. (2018). Gazepath: An eye-tracking analysis tool that accounts for individual differences and data quality. Behavior Research Methods, 50, 834-852.


# Disclaimer
This software has been adapted for current use, but I am not the original author. It is provided "as is", without any warranties. Use at your own risk. The authors are not liable for any issues or damages resulting from its use.



