# Gazepath updated fucntions - July 8th 2024
This repository is a fork that recreates and updates the gazepath functions for parsing gaze data into fixations and saccades, as the original gazepath package is no longer available and does not run on newer versions of R. I have adapted and updated the functions to run directly in R without the need for the package.


# Example of Usage
run the code u the scripts ```helperfunctions.R``` and  ```GazePathHC.R``` to recreate the necessary fucntions


```
dta_4gazepath_stim <- as.data.frame(as.matrix(dta_4gazepath_stim))
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

# Example plots that compare before and after to make sure that parsed data matched the raw data
tmp_stim_by_pp$`140` %>%
  subset(trial_no == 6) %>%
  ggplot(aes(as.numeric(gaze_x_cor_pix), as.numeric(gaze_y_cor_pix))) +
  geom_line() +
  geom_line(aes(y = as.numeric(gaze_y_cor_pix)))

tmp_stim_by_pp$`140` %>%
  subset(trial_no == 6) %>%
  ggplot(aes(as.numeric(timerezero), as.numeric(gaze_x_cor_pix))) +
  geom_line() +
  geom_line(aes(y = as.numeric(gaze_y_cor_pix)))

plot.gazepath(gazepath_results$`140`, trial_index = 6)

```
To Do
Organise the code and fix typos
Split helper functions into individual scripts
Recreate the Shiny app

If you use this, please cite the original authors:
van Renswoude, D. R., Raijmakers, M. E., Koornneef, A., Johnson, S. P., Hunnius, S., & Visser, I. (2018). Gazepath: An eye-tracking analysis tool that accounts for individual differences and data quality. Behavior Research Methods, 50, 834-852.

