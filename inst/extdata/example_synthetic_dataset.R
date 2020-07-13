# Example pipeline
# First, the example dataset stored in synthetic_dataset_alpha.dat is analyzed and used to train a Bayesian classifier.
# Second, the example dataset stored in synthetic_dataset_beta.dat is analyzed and classified.

# First of all, load the MiRNA-QC-and-Diagnosis package:
library(MiRNAQCD)

###########################
### Classifier training ###
###########################

### Loading dataset alpha
# A simple "read.table" does the job:
datasetAlpha <- read.table(file="dataset_alpha.dat", header=TRUE)

### Preprocessing
# The miRNA_expressionPreprocessing function computes sample mean and standard deviation of multiplets.
preprocessedDatasetAlpha <- miRNA_expressionPreprocessing(datasetAlpha, multipletSize=3)

### Outlier removal
# The miRNA_assessQualityThreshold function assesses the critical standard deviation for outlier removal.
qualityThresholdAlpha <- miRNA_assessQualityThreshold(preprocessedDatasetAlpha, significanceLevel=0.05, saveOutputFile=TRUE, outputFileName="qualityThreshold_datasetAlpha.dat")
# The miRNA_removeOutliers function removes any outlier from a dataset according to the assessed critical values.
purifiedDatasetAlpha <- miRNA_removeOutliers(preprocessedDatasetAlpha, qualityThresholdAlpha)

### Feature analysis
# The Target and Versus sets first have to be defined.
Target <- c("A")
Versus <- c("B", "C")
# The miRNA_classifierSetup function, without any feature list, runs in 'Analysis mode' and carries out the analysis of all features.
statisticsAlpha <- miRNA_classifierSetup(purifiedDatasetAlpha, inputTargetList=Target, inputVersusList=Versus, saveOutputFile=TRUE, outputFileBasename="mirnaAnalysis_datasetAlpha")
# The same as the last call, but using "FZ" as normalizer.
statisticsAlphaNorm <- miRNA_classifierSetup(purifiedDatasetAlpha, inputTargetList=Target, inputVersusList=Versus, inputMiRNAList=c("FZ"), saveOutputFile=TRUE, outputFileBasename="mirnaAnalysisNorm_datasetAlpha")
# The analysis reveals a possible set of features/coefficients to be used for the classification:
mirnaToUse <- c("FX", "FZ")
coefficientsToUse <- c(1.0, -1.0)

### Training of Bayesian classifier
# The miRNA_classifierSetup function is run in 'Training mode'.
thresholdValues <- miRNA_classifierSetup(purifiedDatasetAlpha, inputTargetList=Target, inputVersusList=Versus, inputMiRNAList=mirnaToUse, coeffList=coefficientsToUse, saveOutputFile=TRUE, outputFileBasename="threshold_datasetAlpha")

###################################
### Classification of a dataset ###
###################################

### Loading dataset beta
datasetBeta <- read.table(file="dataset_beta.dat", header=TRUE)

### Preprocessing
preprocessedDatasetBeta <- miRNA_expressionPreprocessing(datasetBeta, multipletSize=3)

# Critical sigma values are loaded from the previously generated file
qualityThresholdValues <- miRNA_loadQualityThreshold("qualityThreshold_datasetAlpha.dat")
# Outlier removal
purifiedDatasetBeta <- miRNA_removeOutliers(preprocessedDatasetBeta, qualityThresholdValues)

### Classification
# Threshold values are loaded from the previously generated file
thresholdValues <- miRNA_loadDiagnosticThreshold("threshold_datasetAlpha.txt")
# Classifier parameters are defined in the same way as for training:
mirnaToUse <- c("FX", "FZ")
coefficientsToUse <- c(1.0, -1.0)
# Classification is carried out by the miRNA_diagnosis function:
diagnosedBeta <- miRNA_diagnosis(purifiedDatasetBeta, inputMiRNAList=mirnaToUse, coeffList=coefficientsToUse, inputThreshold=thresholdValues, saveOutputFile=TRUE, outputFileBasename="diagnosis_datasetBeta")
