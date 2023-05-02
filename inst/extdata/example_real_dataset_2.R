########################################
### Example pipeline, Real dataset 2 ###
########################################

## This pipeline reproduces the results of the work
# S. Detassis, V. del Vescovo, M. Grasso, S. Masella, C. Cantaloni, L. Cima, A. Cavazza, P. Graziano, G. Rossi, M. Barbareschi, L. Ricci and M. A. Denti, "miR375-3p Distinguishes Low-Grade Neuroendocrine From Non-neuroendocrine Lung Tumors in FFPE Samples", Frontiers in Molecular Biosciences 7:86 (2020) doi:10.3389/fmolb.2020.00086

## Data concerning this work have been submitted to the GEO database and a link for direct download will be provided as soon as it is available.

## A copy of the data is stored in two files within the same directory of the present script
#	real_dataset_2_training.dat	Training set
#	real_dataset_2_testing.dat	Testing set

#########################################
### Load data and prepare data frames ###
#########################################

# First of all, load the MiRNA-QC-and-Diagnosis package:
library(MiRNAQCD)

# A simple "read.table" does the job
rawTrainingSet = read.table(file="real_dataset_2_training.dat", fileEncoding="UTF-8", header=FALSE)
rawTestingSet = read.table(file="real_dataset_2_testing.dat", fileEncoding="UTF-8", header=FALSE)
# In this case, no headers are loaded from file, so column names must be properly set
names(rawTrainingSet) <- c("miRNA", "Subject", "Value", "Class")
names(rawTestingSet) <- c("miRNA", "Subject", "Value", "Class")

# For the sake of keeping this directory clean, a subdirectory is created to store results
dir.create(file.path(getwd(), "Results_real_dataset_2"), showWarnings = FALSE)

#################################################
### Analysis using MiRNAQCD package functions ###
#################################################

### Preprocessing
# The miRNA_expressionPreprocessing function computes sample mean and standard deviation of multiplets.
preprocTrainingSet <- miRNA_expressionPreprocessing(rawTrainingSet, multipletSize=3)
preprocTestingSet <- miRNA_expressionPreprocessing(rawTestingSet, multipletSize=3)

### Outlier removal
# The miRNA_assessQualityThreshold function assesses the critical standard deviation for outlier removal.
qualityThresholdTrainingSet <- miRNA_assessQualityThreshold(preprocTrainingSet, significanceLevel=0.05)
qualityThresholdTestingSet <- miRNA_assessQualityThreshold(preprocTestingSet, significanceLevel=0.05)
# The miRNA_removeOutliers function removes any outlier from a dataset according to the assessed critical values.
cleanedTrainingSet <- miRNA_removeOutliers(preprocTrainingSet, qualityThresholdTrainingSet)
cleanedTestingSet <- miRNA_removeOutliers(preprocTestingSet, qualityThresholdTestingSet)

### Feature analysis
# The Target and Versus sets first have to be defined.
Target <- c("AT","TC")
Versus <- c("AD","SQC")
# The miRNA_classifierSetup function, without any miRNA list, runs in 'Analysis mode' and carries out the analysis of all miRNAs.
# Training dataset, normalizer U6
outputAnalysis <- miRNA_classifierSetup(cleanedTrainingSet, inputTargetList=Target, inputVersusList=Versus, inputMiRNAList="U6", saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_2/miRNAanalysis_training_dataset_norm_U6")

###############################################
### Actual training / testing of classifier ###
###############################################

### Bayesian classifier
# Classifier relies on the linear combination -miR375+U6, as in the reference (see beginning of this file).
outputTraining <- miRNA_classifierSetup(cleanedTrainingSet, inputTargetList=Target, inputVersusList=Versus, inputMiRNAList=c("miR375","U6"), coeffList=c(-1.0, 1.0), histogramParameters="-6.0_14.0_1.0", saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_2/classifier_training_dataset_-miR375_U6", scorePlotAscending=FALSE, scorePlotParameters="-5.0_15.0_5", colorComplementFlag=TRUE)

### Classifier testing
# Test Bayesian classifier using the classifier -miR375+U6
classifiedDataset <- miRNA_diagnosis(cleanedTestingSet, inputMiRNAList=c("miR375","U6"), coeffList=c(-1.0, 1.0), inputThreshold=outputTraining, inputTargetList=Target, inputVersusList=Versus, saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_2/diagnosis_testing_dataset_-miR375_U6", scorePlotAscending=FALSE, scorePlotParameters="-10.0_10.0_5", colorComplementFlag=TRUE)
# Repeat the test by keeping outliers (to reproduce the plots present in the reference)
classifiedDataset <- miRNA_diagnosis(preprocTestingSet, inputMiRNAList=c("miR375","U6"), coeffList=c(-1.0, 1.0), inputThreshold=outputTraining, inputTargetList=Target, inputVersusList=Versus, saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_2/diagnosis_testing_dataset_with_outliers_-miR375_U6", scorePlotAscending=FALSE, scorePlotParameters="-10.0_10.0_5", colorComplementFlag=TRUE)
