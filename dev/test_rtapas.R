library(rtapas)
library(neurobase)
library(oro.nifti)
source('https://raw.githubusercontent.com/avalcarcel9/aliviateR/master/R/dsc.R')
library(dplyr)
library(tibble)


## -----------------------------------------------------------------------------
# Make a list of the gold standard masks
train_gold_standard_masks = list(gs1  = rtapas::gs1,
                                 gs2  = rtapas::gs2,
                                 gs3  = rtapas::gs3,
                                 gs4  = rtapas::gs4,
                                 gs5  = rtapas::gs5,
                                 gs6  = rtapas::gs6,
                                 gs7  = rtapas::gs7,
                                 gs8  = rtapas::gs8,
                                 gs9  = rtapas::gs9,
                                 gs10 = rtapas::gs10)

# Convert the gold standard masks to nifti objects
train_gold_standard_masks = lapply(train_gold_standard_masks, oro.nifti::nifti)

## -----------------------------------------------------------------------------
# Make a list of the training probability maps
train_probability_maps = list(pmap1 = pmap1, 
                              pmap2 = pmap2, 
                              pmap3 = pmap3, 
                              pmap4 = pmap4, 
                              pmap5 = pmap5, 
                              pmap6 = pmap6, 
                              pmap7 = pmap7, 
                              pmap8 = pmap8, 
                              pmap9 = pmap9, 
                              pmap10 = pmap10)

# Convert the probability maps to nifti objects
train_probability_maps = lapply(train_probability_maps, oro.nifti::nifti)

## -----------------------------------------------------------------------------
# Make a list of the brain masks
train_brain_masks = list(brain_mask1 = brain_mask, 
                         brain_mask2 = brain_mask, 
                         brain_mask3 = brain_mask, 
                         brain_mask4 = brain_mask, 
                         brain_mask5 = brain_mask, 
                         brain_mask6 = brain_mask, 
                         brain_mask7 = brain_mask, 
                         brain_mask8 = brain_mask, 
                         brain_mask9 = brain_mask, 
                         brain_mask10 = brain_mask)

# Convert the brain masks to nifti objects
train_brain_masks = lapply(train_brain_masks, oro.nifti::nifti)


## -----------------------------------------------------------------------------
train_ids = paste0('subject_', 1:length(train_gold_standard_masks))


## -----------------------------------------------------------------------------
# Initialize an empty list to store results created on 1 core
train_data1 = list()

# Run tapas_data function
for(i in 1:length(train_probability_maps)){
  train_data1[[i]] = tapas_data(thresholds = seq(from = 0, to = 1, by = 0.01),
                         pmap = train_probability_maps[[i]],
                         gold_standard = train_gold_standard_masks[[i]],
                         mask = train_brain_masks[[i]],
                         k = 0,
                         subject_id = train_ids[[i]],
                         verbose = TRUE)
}


## -----------------------------------------------------------------------------
# Bind the single core data across list objects (subjects)
train_data1 = dplyr::bind_rows(train_data1)

## -----------------------------------------------------------------------------
tapas_model = tapas_train(data = train_data1, 
                          dsc_cutoff = 0.03, 
                          verbose = TRUE)


## -----------------------------------------------------------------------------
# The TAPAS GAM model
summary(tapas_model$tapas_model)
# The threshold that optimizes group-level DSC
tapas_model$group_threshold
# The lower and upper bound clamps to avoid extrapolation
tapas_model$clamp_data
# The training data for the TAPAS `mgcv::gam` function
tapas_model$train_data


## -----------------------------------------------------------------------------
#make_scattergram(tapas_model = tapas_model)


## -----------------------------------------------------------------------------
# Testing gold standard masks
test_gold_standard_masks = list(gs11 = gs11,
                                gs12 = gs12,
                                gs13 = gs13,
                                gs14 = gs14,
                                gs15 = gs15)

# Make array objects niftis
test_gold_standard_masks = lapply(test_gold_standard_masks, oro.nifti::nifti)

## -----------------------------------------------------------------------------
# Obtain the test subject probability maps
test_probability_maps = list(pmap11 = pmap11, 
                             pmap12 = pmap12,
                             pmap13 = pmap13,
                             pmap14 = pmap14,
                             pmap15 = pmap15)

# Make array objects niftis
test_probability_maps = lapply(test_probability_maps, oro.nifti::nifti)

## -----------------------------------------------------------------------------
# Create a list of testing brain masks
test_brain_masks = list(brain_mask11 = brain_mask, 
                        brain_mask12 = brain_mask,
                        brain_mask13 = brain_mask,
                        brain_mask14 = brain_mask,
                        brain_mask15 = brain_mask)

# Make array objects niftis
test_brain_masks = lapply(test_brain_masks, oro.nifti::nifti)


## -----------------------------------------------------------------------------
test_ids = paste0('subject_', (10 + 1:length(test_gold_standard_masks)))


## -----------------------------------------------------------------------------
# Initialize an empty list to store results created on 1 core
test_data1 = list()

# Run tapas_predict function
for(i in 1:length(test_probability_maps)){
  test_data1[[i]] = tapas_predict(pmap = test_probability_maps[[i]],
                                  model = tapas_model,
                                  clamp = TRUE,
                                  k = 0,
                                  verbose = TRUE)
}


## -----------------------------------------------------------------------------
# Calculate DSC in each mask
dsc = tibble::tibble(
  tapas_dsc = c(dsc(test_gold_standard_masks[[1]], test_data1[[1]]$tapas_binary_mask),
                dsc(test_gold_standard_masks[[2]], test_data1[[2]]$tapas_binary_mask),
                dsc(test_gold_standard_masks[[3]], test_data1[[3]]$tapas_binary_mask),
                dsc(test_gold_standard_masks[[4]], test_data1[[4]]$tapas_binary_mask),
                dsc(test_gold_standard_masks[[5]], test_data1[[5]]$tapas_binary_mask)),
  group_dsc = c(dsc(test_gold_standard_masks[[1]], test_data1[[1]]$group_binary_mask),
                dsc(test_gold_standard_masks[[2]], test_data1[[2]]$group_binary_mask),
                dsc(test_gold_standard_masks[[3]], test_data1[[3]]$group_binary_mask),
                dsc(test_gold_standard_masks[[4]], test_data1[[4]]$group_binary_mask),
                dsc(test_gold_standard_masks[[5]], test_data1[[5]]$group_binary_mask)))
# Print DSC
dsc
