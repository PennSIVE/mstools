## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ---- eval = FALSE------------------------------------------------------------
## source("https://neuroconductor.org/neurocLite.R")
## neuro_install("mimosa")


## ---- eval = FALSE------------------------------------------------------------
## devtools::install_github('avalcarcel9/mimosa')


## ---- warning = FALSE, message = FALSE----------------------------------------
library(neurobase)
library(mimosa)
library(dplyr)
library(oasis)
library(fslr)


## ---- echo = FALSE------------------------------------------------------------
info = Sys.info()
user = info[["user"]]
if (grepl("muschel", user )) {
  data_dir = "~/Desktop/"
  train_data_dir = data_dir
  test_data_dir = data_dir
} else {
  data_dir = "/Users/alval/Documents/ISBI_Challenge_2015/"
  train_data_dir = file.path(data_dir, "training")
  test_data_dir = file.path(data_dir, "testdata_website")
}
train_dir = file.path(data_dir, "training01")


## -----------------------------------------------------------------------------
# Note these paths will be to where you 
# have stored the data or to your own data
train_dir = file.path(train_data_dir, "training01")

T1_files = list.files(path = 
                        file.path(train_dir,
                                  "preprocessed"), 
                      pattern = "mprage_pp[.]nii", 
                      full.names = TRUE)
T2_files = list.files(path = 
                        file.path(train_dir,
                                  "preprocessed"), 
                      pattern = "t2_pp[.]nii", 
                      full.names = TRUE)
FLAIR_files = list.files(path = 
                           file.path(train_dir,
                                     "preprocessed"), 
                         pattern = "flair_pp[.]nii", 
                         full.names = TRUE)
PD_files = list.files(path = 
                        file.path(train_dir,
                                  "preprocessed"), 
                      pattern = "pd_pp[.]nii", 
                      full.names = TRUE)
GS_files = list.files(path = 
                        file.path(train_dir,
                                  "masks"), 
                      pattern = "mask2[.]nii", 
                      full.names = TRUE)
filepaths = data.frame(T1 = T1_files, T2 = T2_files, 
                       FLAIR = FLAIR_files, PD = PD_files, GS = GS_files,
                       stringsAsFactors = FALSE)
have_data = nrow(filepaths) > 0
if (have_data) {
  ss = strsplit(nii.stub(filepaths$T1), split = "_")
  filepaths$visit_id = sapply(ss, function(x) x[2])
  filepaths
}


## ---- include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(eval = have_data)


## ---- warning = FALSE---------------------------------------------------------
# The R package neurobase is needed for the readnii function
T1_training01_01 = readnii(filepaths$T1[1])
T2_training01_01 = readnii(filepaths$T2[1])
FLAIR_training01_01 = readnii(filepaths$FLAIR[1])
PD_training01_01 = readnii(filepaths$PD[1])
gold_standard = readnii(filepaths$GS[1])


## -----------------------------------------------------------------------------
sequence = list(FLAIR = FLAIR_training01_01,
                T1 = T1_training01_01,
                T2 = T2_training01_01,
                PD = PD_training01_01)
multi_overlay(sequence,
              z = floor(oro.nifti::nsli(sequence[[1]])/2),
              text = names(sequence),
              text.y = rep(1.4, length(sequence)),
              text.cex = rep(2.5, length(sequence))
)
rm(sequence)


## -----------------------------------------------------------------------------
create_brain_mask = function(...) {
  x = list(...)
  x = check_nifti(x)
  x = lapply(x, function(img) {
    img > 0
  })
  mask = Reduce("|", x)
  mask = datatyper(mask)
  mask
}


## -----------------------------------------------------------------------------
# Create a brain mask
brain_mask = create_brain_mask(
  T1_training01_01, 
  T2_training01_01,
  FLAIR_training01_01,
  PD_training01_01
)

# The mimosa R package is needed to run mimosa_data
mimosa_data = mimosa_data(
  brain_mask = brain_mask, 
  FLAIR = FLAIR_training01_01, 
  T1 = T1_training01_01, 
  T2 = T2_training01_01, 
  PD = PD_training01_01, 
  tissue = FALSE, 
  gold_standard = gold_standard, 
  normalize = 'Z', 
  cand_mask = NULL, 
  slices = NULL, 
  orientation = c("axial", "coronal", "sagittal"), 
  cores = 1, 
  verbose = TRUE)



## -----------------------------------------------------------------------------
is.list(mimosa_data)
names(mimosa_data)
names(mimosa_data$smoothed)
names(mimosa_data$smoothed$smooth_10)
names(mimosa_data$smoothed$smooth_20)
names(mimosa_data$coupling_intercepts)
names(mimosa_data$coupling_slopes)
names(mimosa_data$normalized)


## -----------------------------------------------------------------------------
head(mimosa_data$mimosa_dataframe)

ortho2(mimosa_data$top_voxels)
ortho2(mimosa_data$smoothed$smooth_10$FLAIR_10)
ortho2(mimosa_data$coupling_slopes$FLAIRonT1_slopes)

# Remove mimosa_data from memory
rm(mimosa_data)


## -----------------------------------------------------------------------------
filepaths$brainmask = NA

# The neurobase R package is required to read and write images
for (i in seq(nrow(filepaths))) {
  # Load files
  visit_id = filepaths$visit_id[i]
  fname = file.path(train_dir, 
                    "preprocessed", 
                    paste0("brainmask_",
                           visit_id, ".nii.gz"))
  if (!file.exists(fname)) {
    T1_training = readnii(filepaths$T1[i])
    T2_training = readnii(filepaths$T2[i])
    FLAIR_training = readnii(filepaths$FLAIR[i])
    PD_training = readnii(filepaths$PD[i])
    brain_mask = create_brain_mask(
      T1_training,
      T2_training,
      FLAIR_training,
      PD_training
    )
  # Save brain mask to local working directory
  writenii(brain_mask, 
           filename = fname)
  }
  filepaths$brainmask[i] = fname
}


## -----------------------------------------------------------------------------
mimosa_training = mimosa_training(
  brain_mask = filepaths$brainmask,
  FLAIR = filepaths$FLAIR,
  T1 = filepaths$T1,
  T2 = filepaths$T2,
  PD = filepaths$PD,
  tissue = FALSE, 
  gold_standard = filepaths$GS,
  normalize = 'Z', 
  slices = NULL, 
  orientation = c("axial", "coronal", "sagittal"),
  cores = 1, 
  verbose = TRUE, 
  outdir = NULL, 
  optimal_threshold = seq(0.25, 0.35, 0.01))
 
names(mimosa_training)
mimosa_training$mimosa_fit_model
mimosa_training$estimated_optimal_threshold


## -----------------------------------------------------------------------------
# Initialize an empty list
mimosa_df_list = vector(mode = "list",
  length = nrow(filepaths))
names(mimosa_df_list) = filepaths$visit_id

for (i in seq(nrow(filepaths))) {
  # Load files
  T1_training = readnii(filepaths$T1[i])
  T2_training = readnii(filepaths$T2[i])
  FLAIR_training = readnii(filepaths$FLAIR[i])
  PD_training = readnii(filepaths$PD[i])
  gold_standard = readnii(filepaths$GS[i])
  brain_mask = readnii(filepaths$brainmask[i])
  # Obtain the mimosa predictor data.frame
  
  mimosa_df_list[[i]] = mimosa_data(
    brain_mask = brain_mask, 
    FLAIR = FLAIR_training, 
    T1 = T1_training,
    T2 = T2_training, 
    PD = PD_training, 
    tissue = FALSE, 
    gold_standard = gold_standard, 
    normalize = 'Z', 
    cand_mask = NULL, 
    slices = NULL, 
    orientation = c("axial", "coronal", "sagittal"), 
    cores = 1, 
    verbose = TRUE)$mimosa_dataframe
}
# Turn list into a single data.frame which has all subjects predictor data.frames
mimosa_df = dplyr::bind_rows(mimosa_df_list, .id = "visit_id")

head(mimosa_df)
dim(mimosa_df)


## -----------------------------------------------------------------------------
formula = gold_standard ~ FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + 
  PD_10 * PD + PD_20 * PD + 
  T2_10 * T2 + T2_20 * T2 + 
  T1_10 * T1 + T1_20 * T1 +
  FLAIRonT1_intercepts + FLAIRonT2_intercepts + FLAIRonPD_intercepts +
  T1onT2_intercepts + T1onPD_intercepts + T2onPD_intercepts +
  T1onFLAIR_intercepts + T2onFLAIR_intercepts + PDonFLAIR_intercepts + 
  T2onT1_intercepts + PDonT1_intercepts + PDonT2_intercepts +
  FLAIRonT1_slopes + FLAIRonT2_slopes + FLAIRonPD_slopes +
  T1onT2_slopes + T1onPD_slopes + T2onPD_slopes +
  T1onFLAIR_slopes + T2onFLAIR_slopes + PDonFLAIR_slopes +
  T2onT1_slopes + PDonT1_slopes + PDonT2_slopes


## -----------------------------------------------------------------------------
formula = gold_standard ~ FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + 
  T2_10 * T2 + T2_20 * T2 + 
  T1_10 * T1 + T1_20 * T1 +
  FLAIRonT1_intercepts + FLAIRonT2_intercepts + 
  T1onT2_intercepts + T1onFLAIR_intercepts + 
  T2onFLAIR_intercepts + T2onT1_intercepts +
  FLAIRonT1_slopes + FLAIRonT2_slopes + 
  T1onT2_slopes + T1onFLAIR_slopes + 
  T2onFLAIR_slopes + T2onT1_slopes


## -----------------------------------------------------------------------------
formula = gold_standard ~ FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + 
  PD_10 * PD + PD_20 * PD + 
  T1_10 * T1 + T1_20 * T1 +
  FLAIRonT1_intercepts + FLAIRonPD_intercepts + 
  T1onPD_intercepts + T1onFLAIR_intercepts + 
  PDonFLAIR_intercepts + PDonT1_intercepts +
  FLAIRonT1_slopes + FLAIRonPD_slopes + 
  T1onPD_slopes + T1onFLAIR_slopes + 
  PDonFLAIR_slopes + PDonT1_slopes


## -----------------------------------------------------------------------------
formula = gold_standard ~ FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + 
  T1_10 * T1 + T1_20 * T1 +
  FLAIRonT1_intercepts + T1onFLAIR_intercepts +
  FLAIRonT1_slopes + T1onFLAIR_slopes


## -----------------------------------------------------------------------------
formula = gold_standard ~ FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + 
  PD_10 * PD + PD_20 * PD + 
  T2_10 * T2 + T2_20 * T2 + 
  T1_10 * T1 + T1_20 * T1 +
  FLAIRonT1_intercepts + FLAIRonT2_intercepts + FLAIRonPD_intercepts +
  T1onT2_intercepts + T1onPD_intercepts + T2onPD_intercepts +
  T1onFLAIR_intercepts + T2onFLAIR_intercepts + PDonFLAIR_intercepts + 
  T2onT1_intercepts + PDonT1_intercepts + PDonT2_intercepts +
  FLAIRonT1_slopes + FLAIRonT2_slopes + FLAIRonPD_slopes +
  T1onT2_slopes + T1onPD_slopes + T2onPD_slopes +
  T1onFLAIR_slopes + T2onFLAIR_slopes + PDonFLAIR_slopes +
  T2onT1_slopes + PDonT1_slopes + PDonT2_slopes

mimosa_model = mimosa_fit(mimosa_df, formula = formula)

mimosa_model


## -----------------------------------------------------------------------------
# The neurobase and mimosa R packages are required for this chunk
# Note these paths will be to where you have stored the data or to your own data
test_dir = file.path(test_data_dir, "test01")

T1_files = list.files(
  path = file.path(test_dir, "preprocessed"), 
                      pattern = "mprage_pp[.]nii", 
                      full.names = TRUE)
T2_files = list.files(
    path = file.path(test_dir, "preprocessed"), 
                      pattern = "t2_pp[.]nii", 
                      full.names = TRUE)
FLAIR_files = list.files(
    path = file.path(test_dir, "preprocessed"), 
                         pattern = "flair_pp[.]nii", 
                         full.names = TRUE)
PD_files = list.files(
    path = file.path(test_dir, "preprocessed"), 
                      pattern = "pd_pp[.]nii", 
                      full.names = TRUE)

filepaths = data.frame(T1 = T1_files, T2 = T2_files, 
  FLAIR = FLAIR_files, PD = PD_files, 
  stringsAsFactors = FALSE)
filepaths

# Load first subject into R
T1_testing = readnii(filepaths$T1[1])
T2_testing = readnii(filepaths$T2[1])
FLAIR_testing = readnii(filepaths$FLAIR[1])
PD_testing = readnii(filepaths$PD[1])

# Create a brain mask
# Create a brain mask
brain_mask = create_brain_mask(
  T1_testing, 
  T2_testing,
  FLAIR_testing,
  PD_testing
)

mimosa_testdata = mimosa_data(
  brain_mask = brain_mask, 
  FLAIR = FLAIR_training, 
  T1 = T1_training,
  T2 = T2_training, 
  PD = PD_training, 
  tissue = FALSE, 
  gold_standard = NULL, 
  normalize = 'Z', 
  cand_mask = NULL, 
  slices = NULL, 
  orientation = c("axial", "coronal", "sagittal"), 
  cores = 1, 
  verbose = TRUE)

mimosa_testdata_df = mimosa_testdata$mimosa_dataframe
mimosa_candidate_mask = mimosa_testdata$top_voxels

rm(T1_files, T2_files, FLAIR_files, PD_files,
  mimosa_testdata)


## -----------------------------------------------------------------------------
# The R package fslr is required to smooth the probability map
predictions = predict(mimosa_model,
                      newdata = mimosa_testdata_df,
                      type = 'response')
probability_map = niftiarr(brain_mask, 0)
probability_map[mimosa_candidate_mask == 1] = predictions

probability_map = fslsmooth(probability_map, 
                            sigma = 1.25,
                            mask = brain_mask, 
                            retimg = TRUE,
                            smooth_mask = TRUE)


## -----------------------------------------------------------------------------
ortho2(probability_map)


## -----------------------------------------------------------------------------
threshold = mimosa_training$estimated_optimal_threshold
segmentation_mask = probability_map > threshold

rm(probability_map)


## -----------------------------------------------------------------------------
ortho2(segmentation_mask)

# The R package scales is needed for a few of the next commands
double_ortho(FLAIR_testing, segmentation_mask, col.y = 'red')
ortho2(FLAIR_testing, segmentation_mask, 
  col.y = "#FF000080")

rm(segmentation_mask)


## -----------------------------------------------------------------------------
train2_dir = file.path(train_data_dir, "training02")

# Read in images
T1_files = list.files(
  path = 
    file.path(train2_dir,
              "preprocessed"), 
  pattern = "mprage_pp[.]nii", 
  full.names = TRUE)
T2_files = list.files(
  path = 
    file.path(train2_dir,
              "preprocessed"), 
  pattern = "t2_pp[.]nii", 
  full.names = TRUE)
FLAIR_files = list.files(
  path = 
    file.path(train2_dir,
              "preprocessed"), 
  pattern = "flair_pp[.]nii", 
  full.names = TRUE)
PD_files = list.files(
  path = 
    file.path(train2_dir,
              "preprocessed"), 
  pattern = "pd_pp[.]nii", 
  full.names = TRUE)
GS_files = list.files(
  path = 
    file.path(train2_dir,
              "masks"), 
  pattern = "mask2[.]nii", 
  full.names = TRUE)
filepaths = data.frame(T1 = T1_files, T2 = T2_files, FLAIR = FLAIR_files, PD = PD_files, GS = GS_files, stringsAsFactors = FALSE)
ss = strsplit(nii.stub(filepaths$T1), split = "_")
filepaths$visit_id = sapply(ss, function(x) x[2])
filepaths

T1 = filepaths$T1[1]
T2 = filepaths$T2[1]
FLAIR = filepaths$FLAIR[1]
PD = filepaths$PD[1]
gold_standard = filepaths$GS[1]

# Create a brain mask
brain_mask = create_brain_mask(
  T1, 
  T2,
  FLAIR,
  PD
)

# Obtain predictor matrix
training02_01_data = mimosa_data(
  brain_mask = brain_mask, 
  FLAIR = FLAIR, 
  T1 = T1,
  T2 = T2, 
  PD = PD, 
  tissue = FALSE, 
  gold_standard = gold_standard, 
  normalize = 'Z', 
  cand_mask = NULL, 
  slices = NULL, 
  orientation = c("axial", "coronal", "sagittal"), 
  cores = 1, 
  verbose = TRUE)

# Create predictions based on trained model
predictions = predict(mimosa_model,
                      newdata = training02_01_data$mimosa_dataframe,
                      type = 'response')

# Create a probability map
probability_map = niftiarr(brain_mask, 0)
probability_map[training02_01_data$top_voxels == 1] = predictions

probability_map = fslsmooth(probability_map, 
                            sigma = 1.25,
                            mask = brain_mask, 
                            retimg = TRUE,
                            smooth_mask = TRUE)

# Threshold probability map to create segmentation mask
threshold = mimosa_training$estimated_optimal_threshold
segmentation_mask = probability_map > threshold

gold_standard = readnii(gold_standard)
# Generate summary measures for performance
count_stats(gold_standard = gold_standard, 
            predicted_segmentation = segmentation_mask, 
            k = 27, 
            percent_overlap = 0.2, 
            verbose = TRUE)

