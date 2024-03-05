# Set working directory
setwd('/Users/jamesboot/Documents/9.Genome Centre Files/Platelet Analysis/Agonist_Only_Analysis/')
base.dir <- '/Users/jamesboot/Documents/9.Genome Centre Files/Platelet Analysis'

# Load packages
library('stringr')
library('ggplot2')
library('dplyr')
library('limma')
library('reshape2')
library('tidyverse')
library('corrplot')
library('cluster') 
library('factoextra')
library('flowCore')
library('FlowSOM')
library('data.table')
library('Spectre')

# Set random seed
set.seed(123)

# Make list of all files
raw.dir <- '/Users/jamesboot/Documents/9.Genome Centre Files/Platelet Analysis/Raw_data/wholeblood'
inFiles <- list.files(path = raw.dir,
                      pattern = '*.fcs',
                      recursive = T,
                      full.names = T)

# Load data using function
all_fcs_files <- lapply(inFiles, function(i){
  read.FCS(i, truncate_max_range = F, which.lines = 200)
})

# Name each element in the loaded file list after the file name 
names(all_fcs_files) <- inFiles

# Check number of cols in each dataframe of list 
ncols_check <- unlist(lapply(all_fcs_files, function(i) {
  ncol(i)
}))
print(paste('How many columns are there across dataframes?', unique(ncols_check)))

# Refine the colnames of all_fcs_files
# Gather all the unique colnames 
cols <- unique(unlist(lapply(all_fcs_files, function(i){
  colnames(i)
})))

# Make a key-values pair list
recol <- c(
  "APC-Cy7-LAMP1-A"='CD41',
  "Alexa Fluor 647-A"='GPVI',
  "BB700-A"='CD62P',
  "BUV395-A"='CD61',
  "BUV737-A"='CD29',
  "BV510-A"='CD9',
  "BV605-A"='CD36',
  "BV650-A"='CD42b',
  "BV711-A"='CLEC-2',
  "BV785-A"='CD31',
  "FITC-A"='PAC-1',
  "PE-A"='CD107a',
  "PE-Cy7-A"='CD63',
  "Syto-13-A"='SYTO',
  "eFluor 450-A"='CD42a',
  "APC-Cy7-A"='CD41',
  "BUV661-A"='CD29',
  'FSC-A'='FSC-A',
  'SSC-A'='SSC-A'
)

# Go through all dataframes in the list and subset and rename
for (x in 1:length(all_fcs_files)) {
  
  # Get the column numbers that we need to subset
  indices <- match(names(recol), colnames(all_fcs_files[[x]]))
  indices <- indices[!is.na(indices)]
  recol2 <- colnames(all_fcs_files[[x]])[indices]
  
  # Select columns that we need and rename
  all_fcs_files[[x]] <- all_fcs_files[[x]][, recol2]
  colnames(all_fcs_files[[x]]) <- recol[match(colnames(all_fcs_files[[x]]),names(recol))]
  
  # Set order of cols
  all_fcs_files[[x]] <- all_fcs_files[[x]][, c("GPVI","CD62P","CD61","CD29","CD9","CD36","CD42b","CLEC-2","CD31","PAC-1",
                                               "CD107a","CD63","CD42a","CD41",'FSC-A','SSC-A')]
  
}

# Check colnames are all matching before the rbind later
colnames_check <- lapply(all_fcs_files, function(i){
  colnames(i)
})

colnames_check <- do.call(rbind.data.frame, colnames_check)
colnames(colnames_check) <- c(1:ncol(colnames_check))

for (x in 1:ncol(colnames_check)) {
  
  message(paste('Do colnames all match for column', x, '?'))
  print(length(unique(c(colnames_check[,x])))==1)
  
}

# Make all list items dataframes
all_fcs_dfs <- list()
for (x in 1:length(all_fcs_files)) {
  all_fcs_dfs[[x]] <- as.data.frame(all_fcs_files[[x]]@exprs)
}
names(all_fcs_dfs) <- names(all_fcs_files)

# Merge all dataframes in list using rbindlist
fcs_dat <- rbindlist(all_fcs_dfs,
                     use.names = T,
                     fill = F,
                     idcol = T)

# Load meta data file
meta_file <- read.csv('meta.csv', header = T)

# Make meta data dataframe
meta_data <- data.frame(cell = c(1:nrow(fcs_dat)),
                        Sample = NA,
                        Treatment = NA,
                        Batch = NA,
                        Donor = NA)

# Populate meta data dataframe
for (x in 1:nrow(meta_file)) {
  meta_data$Sample[fcs_dat$.id == meta_file$File[x]] <- meta_file$Sample[x]
  meta_data$Treatment[fcs_dat$.id == meta_file$File[x]] <- meta_file$Treatment[x]
  meta_data$Batch[fcs_dat$.id == meta_file$File[x]] <- meta_file$Batch[x]
  meta_data$Donor[fcs_dat$.id == meta_file$File[x]] <- meta_file$Donor[x]
}

# Remove .id column from fcs_dat
fcs_dat$.id <- NULL

# Next transform the values
fcs_dat_trans <- do.logicle(fcs_dat,
                            colnames(fcs_dat))

# Subset down to transformed data
fcs_dat_trans <- as.data.frame(fcs_dat_trans[,c(17:32)])

# Before proceeding do some data checks
# Check what the minimum value is in each column
fcs_dat_trans %>% summarise_if(is.numeric, min)

# Count the number of negative values in a column
proportion_negatives <- data.frame(row.names = colnames(fcs_dat_trans))
for (x in 1:16) {
  n_negs <- sum(fcs_dat_trans[, x] < 0)
  n_total <- nrow(fcs_dat_trans)
  proportion_neg <- (n_negs / n_total) * 100
  proportion_negatives[x, 1] <- proportion_neg
}
proportion_negatives

# Now data is loaded and pre-processed - start exploratory analysis here ----

# Analyse all agonists together without Vehicle

# Subset the loaded data frame to contain everything except VEH
# Select rows based on treatment meta
fcs_dat_sub <- fcs_dat_trans[!(meta_data$Treatment == 'VEH'), ]
row.names(fcs_dat_sub)

# Subset the meta
meta_data_sub <- meta_data[!(meta_data$Treatment == 'VEH'), ]

# We will need a markers vector
markers <- colnames(fcs_dat_sub)

# Now perform analysis using Sepctre functions
fcs_dat_sub <- as.data.table(fcs_dat_sub)

# Clustering
cyt <- run.flowsom(fcs_dat_sub,
                   use.cols = colnames(fcs_dat_sub))

# tSNE
cyt <- run.tsne(cyt,
                use.cols = colnames(fcs_dat_sub))
  
# Plot clusters
make.colour.plot(
  cyt,
  "tSNE_X",
  "tSNE_Y",
  "FlowSOM_metacluster",
  col.type = 'factor',
  add.label = TRUE,
  dot.size = 0.2,
  title = 'FlowSOM Clusters',
  save.to.disk = F
) +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = 'Cluster'))

ggsave(
  paste0('Ag_Only_Cluster.png'),
  height = 2099,
  width = 2099,
  dpi = 300,
  units = 'px'
)
  
# Heatmap of cluster marker expression
exp <- do.aggregate(cyt,
                    colnames(fcs_dat_sub),
                    by = "FlowSOM_metacluster")
make.pheatmap(
  exp,
  "FlowSOM_metacluster",
  colnames(fcs_dat_sub),
  file.name = paste0('Ag_Only_cluster_heatmap.png')
)
dev.off()
  
# Plot treatment
# Add treatment to cyt
cyt$Treatment <- meta_data_sub$Treatment
make.colour.plot(
  cyt,
  "tSNE_X",
  "tSNE_Y",
  "Treatment",
  col.type = 'factor',
  add.label = TRUE,
  dot.size = 0.2,
  save.to.disk = F
) +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = 'Treatment'))

ggsave(
  paste0('Ag_Only_Treatment.png'),
  plot = last_plot(),
  height = 2099,
  width = 2099,
  dpi = 300,
  units = 'px'
)
  
# Plot batch
# Add batch to cyt
cyt$Batch <- meta_data_sub$Batch
make.colour.plot(
  cyt,
  "tSNE_X",
  "tSNE_Y",
  "Batch",
  col.type = 'factor',
  add.label = F,
  dot.size = 0.2,
  save.to.disk = F
) +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = 'Batch'))

ggsave(
  paste0('Ag_Only_Batch.png'),
  plot = last_plot(),
  height = 2099,
  width = 2099,
  dpi = 300,
  units = 'px'
)
  

  
# Plot donor
# Add donor to cyt
cyt$Donor <- meta_data_sub$Donor
make.colour.plot(
  cyt,
  "tSNE_X",
  "tSNE_Y",
  "Donor",
  col.type = 'factor',
  add.label = F,
  dot.size = 0.2,
  save.to.disk = F
) +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = 'Donor'))

ggsave(
  paste0('Ag_Only_Donor.png'),
  plot = last_plot(),
  height = 2099,
  width = 2099,
  dpi = 300,
  units = 'px'
)
  
# Summary of sample distrbution etc.
cyt$Sample <- meta_data_sub$Sample
# Keep only relevant data
post_meta <-
  cyt[, c('Treatment', 'FlowSOM_metacluster', 'Donor', 'Batch')]

# Summarise using dplyr
summary <- post_meta %>%
  group_by(FlowSOM_metacluster) %>%
  count(Donor) %>%
  mutate(proportion = n / sum(n)) %>%
  mutate(
    Cluster = case_when(
      FlowSOM_metacluster == 1 ~ 'One',
      FlowSOM_metacluster == 2 ~ 'Two',
      FlowSOM_metacluster == 3 ~ 'Three',
      FlowSOM_metacluster == 4 ~ 'Four',
      FlowSOM_metacluster == 5 ~ 'Five',
      FlowSOM_metacluster == 6 ~ 'Six',
      FlowSOM_metacluster == 7 ~ 'Seven',
      FlowSOM_metacluster == 8 ~ 'Eight',
      FlowSOM_metacluster == 9 ~ 'Nine',
      FlowSOM_metacluster == 10 ~ 'Ten',
      FlowSOM_metacluster == 11 ~ 'Eleven',
      FlowSOM_metacluster == 12 ~ 'Twelve',
      FlowSOM_metacluster == 13 ~ 'Thirteen',
      FlowSOM_metacluster == 14 ~ 'Fourteen',
      FlowSOM_metacluster == 15 ~ 'Fifteen',
      FlowSOM_metacluster == 16 ~ 'Sixteen'
    )
  )

summary$Cluster <-
  factor(
    summary$Cluster,
    levels = c(
      'One',
      'Two',
      'Three',
      'Four',
      'Five',
      'Six',
      'Seven',
      'Eight',
      'Nine',
      'Ten',
      'Eleven',
      'Twelve',
      'Thirteen',
      'Fourteen',
      'Fifteen',
      'Sixteen'
    )
  )

# Plot summary stacked bar chart
ggplot(summary, aes(x = Cluster, y = n, fill = Donor)) +
  geom_bar(stat = "identity", position = 'stack') +
  guides(fill = guide_legend(title = "Donor")) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ))

ggsave(
  plot = last_plot(),
  filename =  paste0('Ag_Only_Donor_Distr.png'),
  height = 2099,
  width = 2099,
  dpi = 300,
  units = 'px'
)

# Make a summary of treatment distribution
# Summarise using dplyr
summary <- post_meta %>%
  group_by(FlowSOM_metacluster) %>%
  count(Treatment) %>%
  mutate(proportion = n / sum(n)) %>%
  mutate(
    Cluster = case_when(
      FlowSOM_metacluster == 1 ~ 'One',
      FlowSOM_metacluster == 2 ~ 'Two',
      FlowSOM_metacluster == 3 ~ 'Three',
      FlowSOM_metacluster == 4 ~ 'Four',
      FlowSOM_metacluster == 5 ~ 'Five',
      FlowSOM_metacluster == 6 ~ 'Six',
      FlowSOM_metacluster == 7 ~ 'Seven',
      FlowSOM_metacluster == 8 ~ 'Eight',
      FlowSOM_metacluster == 9 ~ 'Nine',
      FlowSOM_metacluster == 10 ~ 'Ten',
      FlowSOM_metacluster == 11 ~ 'Eleven',
      FlowSOM_metacluster == 12 ~ 'Twelve',
      FlowSOM_metacluster == 13 ~ 'Thirteen',
      FlowSOM_metacluster == 14 ~ 'Fourteen',
      FlowSOM_metacluster == 15 ~ 'Fifteen',
      FlowSOM_metacluster == 16 ~ 'Sixteen'
    )
  )

summary$Cluster <-
  factor(
    summary$Cluster,
    levels = c(
      'One',
      'Two',
      'Three',
      'Four',
      'Five',
      'Six',
      'Seven',
      'Eight',
      'Nine',
      'Ten',
      'Eleven',
      'Twelve',
      'Thirteen',
      'Fourteen',
      'Fifteen',
      'Sixteen'
    )
  )

# Plot summary stacked bar chart
ggplot(summary, aes(x = Cluster, y = n, fill = Treatment)) +
  geom_bar(stat = "identity", position = 'stack') +
  guides(fill = guide_legend(title = "Treatment")) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ))

ggsave(
  plot = last_plot(),
  filename = 'Ag_Only_Treatment_Distr.png',
  height = 2099,
  width = 2099,
  dpi = 300,
  units = 'px'
)

# Now start ML analysis here: just two agonists together ----

# Load package
library(caret)

# Subset down to the selected treatments - ADP plus U4 (weak and strong)
# For both data and metadata
data_new <- fcs_dat_trans[meta_data$Treatment %in% c('ADP', 'U4'), ]
meta_new <- meta_data[meta_data$Treatment %in% c('ADP', 'U4'), ]

# Perform exploratory analysis of these two treatments together

# Clustering - use all columns
cytNew <- run.flowsom(as.data.table(data_new),
                      use.cols = colnames(data_new))

# tSNE - use all columns
cytNew <- run.tsne(cytNew,
                   use.cols = colnames(data_new))

# Visualise results with treatment
# Add treatment to cytNew
cytNew$Treatment <- meta_new$Treatment
make.colour.plot(
  cytNew,
  "tSNE_X",
  "tSNE_Y",
  "Treatment",
  col.type = 'factor',
  add.label = TRUE,
  dot.size = 0.2,
  save.to.disk = F
) +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = 'Treatment'))

# Now create training dataset and unseen validation dataset

# Use 4 donors as training and 2 as unseen validation
# Make the selection at random
trainingDonors <- sample(unique(meta_data$Donor), size = 4, replace = F)

# Add identity column to the dataset
data_new$identity <- as.factor(meta_new$Treatment)

# Subset data 
dataTraining <-
  data_new[(meta_new$Donor %in% trainingDonors), ]
dataUnseen <-
  data_new[(! meta_new$Donor %in% trainingDonors), ]

# Subset meta
metaTraining <-
  meta_new[(meta_new$Donor %in% trainingDonors), ]
metaUnseen <-
  meta_new[(! meta_new$Donor %in% trainingDonors), ]

# Check there are equal numbers of each treatment in training and unseen dataset
sum(dataTraining$identity == 'ADP')
sum(dataTraining$identity == 'U4')
sum(dataUnseen$identity == 'ADP')
sum(dataUnseen$identity == 'U4')

# Create a list of 80% of the rows in the training dataset we can use for training
validation_index <-
  createDataPartition(dataTraining$identity, p = 0.80, list = FALSE)

# Select 20% of the data for initial validation
initialValidation <- dataTraining[-validation_index, ]

# Use the remaining 80% of data to training and testing the models
finalTraining <- dataTraining[validation_index, ]

# Define training parameters
# Run algorithms using 10-fold cross validation repeated 3 times
control <-
  trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 3,
    savePredictions = T,
    classProbs = T,
    summaryFunction = twoClassSummary
  )

# Set metric for performance
metric <- "ROC"

# Train models

# Random Forest
set.seed(7)
fit.rf <-
  train(
    identity ~ .,
    data = finalTraining,
    method = "rf",
    metric = metric,
    trControl = control
  )

# Gradient boosting
set.seed(7)
fit.gbm <-
  train(
    identity ~ .,
    data = finalTraining,
    method = "gbm",
    metric = metric,
    trControl = control
  )

# Boosted generalised linear model
set.seed(7)
fit.glmboost <-
  train(
    identity ~ .,
    data = finalTraining,
    method = "glmboost",
    metric = metric,
    trControl = control
  )

# Support vector machines
set.seed(7)
fit.svmPoly <-
  train(
    identity ~ .,
    data = finalTraining,
    method = "svmPoly",
    metric = metric,
    trControl = control
  )

# K nearest neighbours
set.seed(7)
fit.knn <-
  train(
    identity ~ .,
    data = finalTraining,
    method = "knn",
    metric = metric,
    trControl = control
  )

# Summarize accuracy of models
results <- resamples(
  list(
    rf = fit.rf,
    gbm = fit.gbm,
    glmboost = fit.glmboost,
    svmPoly = fit.svmPoly,
    knn = fit.knn
  )
)

# Compare accuracy of models
dotplot(results)

# Estimate skill of best model on the validation dataset
predictions <- predict(fit.svmPoly, initialValidation)
initialValidResults <- confusionMatrix(predictions, as.factor(initialValidation$identity))

initialValidDF <- data.frame(
  Accuracy = initialValidResults$overall[1],
  Accuracy_Lower = initialValidResults$overall[3],
  Accuracy_Upper = initialValidResults$overall[4],
  Accuracy_Pval = initialValidResults$overall[5],
  row.names = NULL
)

write.csv(initialValidDF, file = 'initialValidationResults.csv')

# Now perform prediction on unseen validation dataset
predictions_unseen <- predict(fit.svmPoly, dataUnseen)
unseenValidResults <- confusionMatrix(predictions_unseen, as.factor(dataUnseen$identity))

unseenValidDF <- data.frame(
  Accuracy = unseenValidResults$overall[1],
  Accuracy_Lower = unseenValidResults$overall[3],
  Accuracy_Upper = unseenValidResults$overall[4],
  Accuracy_Pval = unseenValidResults$overall[5],
  row.names = NULL
)

write.csv(unseenValidDF, file = 'unseenValidationResults.csv')

# Extract variable importance
db.Imp <- varImp(fit.svmPoly, scale = F)
plot(db.Imp)

# Now start ML analysis here: everything together ----

# Subset down to the selected treatments - ADP plus U4 (weak and strong)
# For both data and metadata
data_new <- fcs_dat_sub
meta_new <- meta_data_sub

# Now create training dataset and unseen validation dataset
# Use 4 donors as training and 2 as unseen validation
# Make the selection at random
trainingDonors <- sample(unique(meta_data$Donor), size = 4, replace = F)

# Add identity column to the dataset
data_new$identity <- as.factor(meta_new$Treatment)

# Subset data 
dataTraining <-
  data_new[(meta_new$Donor %in% trainingDonors), ]
dataUnseen <-
  data_new[(! meta_new$Donor %in% trainingDonors), ]

# Subset meta
metaTraining <-
  meta_new[(meta_new$Donor %in% trainingDonors), ]
metaUnseen <-
  meta_new[(! meta_new$Donor %in% trainingDonors), ]

# Check there are equal numbers of each treatment in training and unseen dataset
unique(dataTraining$identity)
sum(dataTraining$identity == 'ADP')
sum(dataTraining$identity == 'U4')
sum(dataUnseen$identity == 'ADP')
sum(dataUnseen$identity == 'U4')

# Create a list of 80% of the rows in the training dataset we can use for training
validation_index <-
  createDataPartition(dataTraining$identity, p = 0.80, list = FALSE)

# Select 20% of the data for initial validation
initialValidation <- dataTraining[-validation_index, ]

# Use the remaining 80% of data to training and testing the models
finalTraining <- dataTraining[validation_index, ]

# Define training parameters
# Run algorithms using 10-fold cross validation repeated 3 times
control <-
  trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 3,
    savePredictions = T,
    classProbs = T,
    summaryFunction = defaultSummary
  )

# Set metric for performance
metric <- "ROC"

# Train models

# Random Forest
set.seed(7)
fit.rf <-
  train(
    identity ~ .,
    data = finalTraining,
    method = "rf",
    metric = metric,
    trControl = control
  )

# Gradient boosting
set.seed(7)
fit.gbm <-
  train(
    identity ~ .,
    data = finalTraining,
    method = "gbm",
    metric = metric,
    trControl = control
  )

# Boosted generalised linear model
set.seed(7)
fit.glmboost <-
  train(
    identity ~ .,
    data = finalTraining,
    method = "glmboost",
    metric = metric,
    trControl = control
  )

# Support vector machines
set.seed(7)
fit.svmPoly <-
  train(
    identity ~ .,
    data = finalTraining,
    method = "svmPoly",
    metric = metric,
    trControl = control
  )

# K nearest neighbours
set.seed(7)
fit.knn <-
  train(
    identity ~ .,
    data = finalTraining,
    method = "knn",
    metric = metric,
    trControl = control
  )

# Summarize accuracy of models
results <- resamples(
  list(
    rf = fit.rf,
    gbm = fit.gbm,
    #glmboost = fit.glmboost,
    svmPoly = fit.svmPoly,
    knn = fit.knn
  )
)

# Compare accuracy of models
dotplot(results)

# Estimate skill of best model on the validation dataset
predictions <- predict(fit.rf, initialValidation)
initialValidResults <- confusionMatrix(predictions, as.factor(initialValidation$identity))

initialValidDF <- data.frame(
  Accuracy = initialValidResults$overall[1],
  Accuracy_Lower = initialValidResults$overall[3],
  Accuracy_Upper = initialValidResults$overall[4],
  Accuracy_Pval = initialValidResults$overall[5],
  row.names = NULL
)

write.csv(initialValidDF, file = 'initialValidationResults.csv')

# Now perform prediction on unseen validation dataset
predictions_unseen <- predict(fit.rf, dataUnseen)
unseenValidResults <- confusionMatrix(predictions_unseen, as.factor(dataUnseen$identity))

unseenValidDF <- data.frame(
  Accuracy = unseenValidResults$overall[1],
  Accuracy_Lower = unseenValidResults$overall[3],
  Accuracy_Upper = unseenValidResults$overall[4],
  Accuracy_Pval = unseenValidResults$overall[5],
  row.names = NULL
)

write.csv(unseenValidDF, file = 'unseenValidationResults.csv')

# Extract variable importance
db.Imp <- varImp(fit.rf, scale = F)
plot(db.Imp)




