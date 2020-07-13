########################################
### Example pipeline, Real dataset 1 ###
########################################

## This pipeline reproduces the results of the work
# L. Ricci, V. Del Vescovo, C. Cantaloni, M. Grasso, M. Barbareschi and M. A. Denti, "Statistical analysis of a Bayesian classifier based on the expression of miRNAs", BMC Bioinformatics 16:287 (2015) doi:10.1186/s12859-015-0715-9

## Data concerning this work are available as a Supplementary Material and can be downloaded at the following URL
# https://static-content.springer.com/esm/art%3A10.1186%2Fs12859-015-0715-9/MediaObjects/12859_2015_715_MOESM1_ESM.txt

## A copy of the same data is stored in three files within the same directory of the present script
#	real_dataset_1_training.dat	Training set
#	real_dataset_1_testing.dat	Testing set
#	real_dataset_1_additional.dat	Additional set

#########################################
### Load data and prepare data frames ###
#########################################

# First of all, load the MiRNA-QC-and-Diagnosis package:
library(MiRNAQCD)

# A simple "read.table" does the job
rawTrainingSet = read.table(file="real_dataset_1_training.dat", fileEncoding="UTF-8", header=FALSE)
rawTestingSet = read.table(file="real_dataset_1_testing.dat", fileEncoding="UTF-8", header=FALSE)
rawAdditionalSet = read.table(file="real_dataset_1_additional.dat", fileEncoding="UTF-8", header=FALSE)
# In this case, no headers are loaded from file, so column names must be properly set
names(rawTrainingSet) <- c("miRNA", "Subject", "Value", "Class")
names(rawTestingSet) <- c("miRNA", "Subject", "Value", "Class")
names(rawAdditionalSet) <- c("miRNA", "Subject", "Value", "Class")

# For the sake of keeping this directory clean, a subdirectory is created to store results
dir.create(file.path(getwd(), "Results_real_dataset_1"), showWarnings = FALSE)

#################################################
### Analysis using MiRNAQCD package functions ###
#################################################

### Preprocessing
# The miRNA_expressionPreprocessing function computes sample mean and standard deviation of multiplets.
preprocessedTrainingSet <- miRNA_expressionPreprocessing(rawTrainingSet, multipletSize=3)
preprocessedTestingSet <- miRNA_expressionPreprocessing(rawTestingSet, multipletSize=3)
preprocessedAdditionalSet <- miRNA_expressionPreprocessing(rawAdditionalSet, multipletSize=3)

### Outlier removal
# The miRNA_assessQualityThreshold function assesses the critical standard deviation for outlier removal.
qualityThresholdTrainingSet <- miRNA_assessQualityThreshold(preprocessedTrainingSet, significanceLevel=0.05, saveOutputFile=FALSE)
qualityThresholdTestingSet <- miRNA_assessQualityThreshold(preprocessedTestingSet, significanceLevel=0.05, saveOutputFile=FALSE)
qualityThresholdAdditionalSet <- miRNA_assessQualityThreshold(preprocessedAdditionalSet, significanceLevel=0.05, saveOutputFile=FALSE)
# The miRNA_removeOutliers function removes any outlier from a dataset according to the assessed critical values.
cleanedTrainingSet <- miRNA_removeOutliers(preprocessedTrainingSet, qualityThresholdTrainingSet)
cleanedTestingSet <- miRNA_removeOutliers(preprocessedTestingSet, qualityThresholdTestingSet)
cleanedAdditionalSet <- miRNA_removeOutliers(preprocessedAdditionalSet, qualityThresholdAdditionalSet)

### Feature analysis
# The Target and Versus sets first have to be defined.
Target <- "ADC"
Versus <- "SQC"
# The miRNA_classifierSetup function, without any miRNA list, runs in 'Analysis mode' and carries out the analysis of all miRNAs.
# Training dataset, no normalizer
outputAnalisysTrainingSet <- miRNA_classifierSetup(cleanedTrainingSet, inputTargetList=Target, inputVersusList=Versus, saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_1/miRNAanalysis_training_dataset")
# Training dataset, normalizer U6
outputAnalisysTrainingSet <- miRNA_classifierSetup(cleanedTrainingSet, inputTargetList=Target, inputVersusList=Versus, inputMiRNAList="U6", saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_1/miRNAanalysis_training_dataset_norm_U6")
# Additional dataset, no normalizer
outputAnalisysAdditionalSet <- miRNA_classifierSetup(cleanedAdditionalSet, inputTargetList=Target, inputVersusList=Versus, saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_1/miRNAanalysis_additional_dataset")
# Additional dataset, normalizer U6
outputAnalisysAdditionalSet_norm <- miRNA_classifierSetup(cleanedAdditionalSet, inputTargetList=Target, inputVersusList=Versus, inputMiRNAList="U6", saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_1/miRNAanalysis_additional_dataset_norm_U6")
# MiRNA analysis is carried out again by selecting either one miRNA + U6 or U6 only.
# This is done to produce histograms as in the reference (see beginning of this file).
trainingSet_dx205 <- cleanedTrainingSet[cleanedTrainingSet$miRNA!="miR21",]
trainingSet_dx21 <- cleanedTrainingSet[cleanedTrainingSet$miRNA!="miR205",]
trainingSet_xU6 <- cleanedTrainingSet[cleanedTrainingSet$miRNA=="U6",]
# Training dataset, dx205
outputAnalysisTrainingSet_dx205 <- miRNA_classifierSetup(trainingSet_dx205, inputTargetList=Target, inputVersusList=Versus, inputMiRNAList="U6", saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_1/miRNAanalysis_training_dataset_dx205", histogramParameters="-10.0_6.0_1.0")
# Training dataset, dx21
outputAnalysisTrainingSet_dx21 <- miRNA_classifierSetup(trainingSet_dx21, inputTargetList=Target, inputVersusList=Versus, inputMiRNAList="U6", saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_1/miRNAanalysis_training_dataset_dx21", histogramParameters="-10.0_0.0_1.0")
# Training dataset, U6 only
outputAnalysisTrainingSet_U6 <- miRNA_classifierSetup(trainingSet_xU6, inputTargetList=Target, inputVersusList=Versus, saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_1/miRNAanalysis_training_dataset_U6", histogramParameters="20.0_31.0_1.0")

###############################################
### Actual training / testing of classifier ###
###############################################

### Bayesian classifier
# "Optimal" classifier used in the reference (see beginning of this file).
outputTraining_yOpt <- miRNA_classifierSetup(cleanedTrainingSet, inputTargetList=Target, inputVersusList=Versus, inputMiRNAList=c("miR205","miR21", "U6"), coeffList=c(1.0, -0.8, -0.2), scorePlotParameters="-2.0_12.0_2.0", histogramParameters="-2.0_12.0_1.0", saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_1/classifier_training_dataset_yOpt")
# Other classifier to compare performance with yOpt
outputTraining_yDV <- miRNA_classifierSetup(cleanedTrainingSet, inputTargetList=Target, inputVersusList=Versus, inputMiRNAList=c("miR205","miR21", "U6"), coeffList=c(1.0, -0.5, -0.5), saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_1/classifier_training_dataset_yDV")
# Using dx205 as classifier
outputTraining_dx205 <- miRNA_classifierSetup(cleanedTrainingSet, inputTargetList=Target, inputVersusList=Versus, inputMiRNAList=c("miR205", "U6"), coeffList=c(1.0, -1.0), saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_1/classifier_training_dataset_dx205", histogramParameters="-12.0_8.0_1.0", scorePlotParameters="-10.0_6.5_2.0")
# Using dx21 as classifier
outputTraining_dx21 <- miRNA_classifierSetup(cleanedTrainingSet, inputTargetList=Target, inputVersusList=Versus, inputMiRNAList=c("miR21", "U6"), coeffList=c(1.0, -1.0), saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_1/classifier_training_dataset_dx21", histogramParameters="-12.0_8.0_1.0", scorePlotParameters="-12.0_0.0_2.0")
# Using -U6 as classifier
outputTraining_u6 <- miRNA_classifierSetup(cleanedTrainingSet, inputTargetList=Target, inputVersusList=Versus, inputMiRNAList="U6", coeffList=-1.0, saveOutputFile=TRUE, outputFileBasename="Results_real_dataset_1/classifier_training_dataset_-U6", histogramParameters="-12.0_8.0_1.0", scorePlotParameters="-30.0_-20.0_2.0")

### Classifier testing
# Test Bayesian classifier using yOpt on cleanedTestingSet
classifiedDataset_yOpt <- miRNA_diagnosis(cleanedTestingSet, inputMiRNAList=c("miR205","miR21", "U6"), coeffList=c(1.0, -0.8, -0.2), outputTraining_yOpt, inputTargetList=Target, inputVersusList=Versus, saveOutputFile=TRUE, scorePlotParameters="-6.0_9.0_2.0", outputFileBasename="Results_real_dataset_1/diagnosis_testing_dataset_yOpt")
