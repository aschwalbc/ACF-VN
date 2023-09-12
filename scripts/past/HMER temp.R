require(hmer)
require(lhs)

# This R-script is a template for you to perform history matching with emulation on a deterministic model. Note that while 
# some sections contain code that is ready to be run, in other sections you will have to adapt the code to your own model.


#############################################################################################################################
######################################## Parameter ranges, targets and design points ########################################
#############################################################################################################################


################################# Create empty lists to store data from the different waves #################################

ems <- list() # ems[[k]] will contain wave-k emulators
wave_data <- list() # wave_data[[k]] will contain the data used to train and validate wave-k emulators
non_imp_pts <- list() # non_imp_pts[[k]] will contain the non-implausible points generated at the end of wave k

# You may also wish to save to file the data from the different waves

################################### Assign the initial parameter ranges to `ranges` ####################################

# For example if you have three parameters, called `param1`, `param2` and `param3`, with ranges (a1,b1), (a2,b2) and (a3,b3), 
# `ranges` would be:
# ranges = list(
#     param1 = c(a1,b1), 
#     param2 = c(a2,b2),
#     param3 = c(a3,b3)
# )


################################################# Define the `targets` list #################################################

# For example if you have two targets, called `target1` and `target2`, with mean value m1 and standard deviation
# sd1 for `target1`, and lower bound lb2 and upper bound ub2 for `target2`, the list `targets` would be: 
# targets = list(
#     target1 = list(val=m1, sigma=sd1), 
#     target2 = c(lb2,ub2)
# )


############################################## Define a latin hypercube design ##############################################

# This can be done through the function `maximinLHS`, which assumes that each parameter is distributed on [0,1]
initial_LHS_training <- lhs::maximinLHS(10 * length(ranges), length(ranges))
initial_LHS_validation <- lhs::maximinLHS(10 * length(ranges), length(ranges))
initial_LHS <- rbind(initial_LHS_training, initial_LHS_validation)
# Adjust each parameter range to be the corrected one (instead of [0,1]) and add columns names to identify the parameters
initial_points <- setNames(data.frame(t(apply(initial_LHS, 1, 
                                              function(x) x*unlist(lapply(ranges, function(x) x[2]-x[1])) + 
                                                  unlist(lapply(ranges, function(x) x[1]))))), names(ranges))


###################### Assign `initial_points` and the correspondent model outputs to `wave_data[[1]]` ######################

# First run the model on the parameter sets in `initial_points`. Put the outputs in a dataframe `initial_results`, having as 
# many rows as `initial_points` and a column for each model output. The columns should be named according to the names in 
# `targets`. Finally bind `initial_points` and `initial_results`:

wave_data[[1]] <- cbind(initial_points, initial_results)


################################# Split `wave_data[[1]]` into training and validation sets ##################################

training <- wave_data[[1]][1:10 * length(ranges),]
validation <- wave_data[[1]][(10 * length(ranges)+1):20 * length(ranges),]



#############################################################################################################################
########################################################### Wave 1 ##########################################################
#############################################################################################################################


############################# Train wave-1 emulators through the function `emulator_from_data` ##############################

ems[[1]] <- emulator_from_data(training, names(targets), ranges)


######################### Validate wave-1 emulators through the function `validation_diagnostics` ###########################

# Plot the three validation tests for each emulator in `ems[[1]]`

vd <- validation_diagnostics(ems[[1]], validation = validation, targets = targets, plt=TRUE)

# If all emulators pass the three tests, go to the next section to generate non-implausible points.
# If, for example, the emulator for `target1` does not pass all three diagnostics, increase its sigma to improve its 
# performance. The code below multiplies the sigma by a factor 1.2: 
# ems[[1]]$`target1` <- ems[[1]]$`target1`$mult_sigma(1.2)
# Then check if the modified emulator passes the diagnostics:
# vd <- validation_diagnostics(ems[[1]][[i]], validation = validation, targets = targets, plt=TRUE)
# If it does, you are done, otherwise you can further increase the sigma, till all three diagnostic tests are successful. 
# Note that you can also choose not to use an emulator at any given wave, if its performance is particularly poor. 
# The following code removes the emulator for `target1` from the list of emulators in the first wave:
# ems[[1]]$`target1` <- NULL

# You can also automate the diagnostic process, with the code below, which does the following:
# 1. Check if there are misclassifications for each emulator (middle column diagnostic in plot above) with the function 
# classification_diag;
# 2. In case of misclassifications, increase the sigma (say by 10%) and go to step 1.
# 3. Once misclassifications have been eliminated, we check how many validation points fail the first of the diagnostics 
# (left column in plot above) with the function comparison_diag, and discard an emulator if it produces too many failures. 
# The code below implements this, removing emulators for which more than 10% of validation points do not pass the first 
# diagnostic.
                                                                
for (j in 1:length(ems[[1]])) {
      misclass <- nrow(classification_diag(ems[[1]][[j]], targets, validation, plt = FALSE))
      while(misclass > 0) {
        ems[[1]][[j]] <- ems[[1]][[j]]$mult_sigma(1.1)
        misclass <- nrow(classification_diag(ems[[1]][[j]], targets, validation, plt = FALSE))
      }
}
bad.ems <- c()
for (j in 1:length(ems[[1]])) {
          bad.model <- nrow(comparison_diag(ems[[1]][[j]], targets, validation, plt = FALSE))
          if (bad.model > floor(nrow(validation)/10)) {
            bad.ems <- c(bad.ems, j)
   }
}
ems[[1]] <- ems[[1]][!seq_along(ems[[1]]) %in% bad.ems]
                                                                
                                                                
######################## Generate non-implausible points through the function `generate_new_runs` ###########################

# The code below generates points that are deemed non-implausible by all wave-1 emulators and assigns the obtained points to
# `non_imp_pts[[1]]`. If in the validation step you decided to discard one or more emulators, then replace `ems[[1]]` with 
# the correct sublist of emulators. We suggest generating 20 points for each input parameter varied, however this can be 
# changed to create fewer or more points.

non_imp_pts[[1]] <- generate_new_design(ems[[1]],  20 * length(ranges), targets, verbose=TRUE)



#############################################################################################################################
###################################################### Subsequent waves #####################################################
#############################################################################################################################

# The code below can be used for all waves after wave 1. Start with k equal to 2 to perform the second wave, then 
# change k to 3 to perform the third wave, and so on.

k <- 2

                                                                

#################### Assign `non_imp_pts[[k-1]]` and the correspondent model outputs to `wave_data[[k]]` ####################

# First run the model on the parameter sets in `non_imp_pts[[k-1]]`. Put the outputs in a dataframe `results_k`, having as 
# many rows as `non_imp_pts[[k-1]]` and a column for each model output. The columns should be named according to the names in 
# `targets`. Finally bind `non_imp_pts[[k-1]]` and `results_k`:

wave_data[[k]] <- cbind(non_imp_pts[[k-1]], results_k)


################################## Split `wave_data[[k]]` into training and validation sets #################################

t_sample <- sample(1:nrow(wave_data[[k]]), round(length(wave_data[[k]][,1])/2))
training <- wave_data[[k]][t_sample,]
validation <- wave_data[[k]][-t_sample,]


############################# Train wave-k emulators through the function `emulator_from_data` ##############################

# From wave 2 on, set `check.ranges=TRUE` to ensure that the new emulators are trained only on the non-implausible space 
# identified in the previous wave.
                                                                
ems[[k]] <- emulator_from_data(training, names(targets), ranges, check.ranges=TRUE)


######################### Validate wave-k emulators through the function `validation_diagnostics` ###########################

# Plot the three validation tests for each of the trained emulators

vd <- validation_diagnostics(ems[[k]], validation = validation, targets = targets, plt=TRUE)

# If all emulators pass the three tests, go to the next section (point proposal). Otherwise try to improve them by increasing 
# their sigmas. As already mentioned above, you can also choose not to use an emulator at any given wave, if its performance 
# is particularly poor. 
                                                                
# You can also automate the diagnostic check, using the code below
                                                                
for (j in 1:length(ems[[k]])) {
      misclass <- nrow(classification_diag(ems[[k]][[j]], targets, validation, plt = FALSE))
      while(misclass > 0) {
        ems[[k]][[j]] <- ems[[k]][[j]]$mult_sigma(1.1)
        misclass <- nrow(classification_diag(ems[[k]][[j]], targets, validation, plt = FALSE))
      }
}
bad.ems <- c()
for (j in 1:length(ems[[k]])) {
          bad.model <- nrow(comparison_diag(ems[[k]][[j]], targets, validation, plt = FALSE))
          if (bad.model > floor(nrow(validation)/10)) {
            bad.ems <- c(bad.ems, j)
   }
}
ems[[k]] <- ems[[k]][!seq_along(ems[[k]]) %in% bad.ems]


######################## Generate non-implausible points through the function `generate_new_runs` ###########################

# The code below generates points, evaluating their implausibility using all emulators generated in all prior waves. If in 
# the validation step you decided to discard one or more emulators, then replace `ems[[k]]` with the correct list of 
# emulators you want to use. When using the `generate_new_runs` function, make sure to pass all emulators trained so far.
# Note that it is important to put the last-wave emulators first, since the `generate_new_design` picks the parameter ranges
# from the first emulator in the list.

non_imp_pts[[k]] <- generate_new_design(c(ems[[k]],ems[[k-1]],...,ems[[1]]),  20 * length(ranges), targets, verbose=TRUE)



#############################################################################################################################
###################################### Comparing non-implausible points across waves ########################################
#############################################################################################################################


############## Show the distribution of the non-implausible space for waves of interest through `wave_points` ###############

# For example, if you want to visualise the distribution of non-implausible points at the end of wave 1, wave 3 and wave 6, 
# you would use the following code:
# wave_points(list(non_imp_pts[[1]], non_imp_pts[[3]], non_imp_pts[[6]]), input_names = names(ranges))


################### Show how non-implausible points perform against the targets using `simulator_plot` #####################

# For example, if you want to compare how non-implausible points from wave 1, wave 3 ad wave 6 perform against the targets,
# you would use the following code:
# simulator_plot(list(wave_data[[1]], wave_data[[3]], wave_data[[6]]), input_names = names(ranges))

















