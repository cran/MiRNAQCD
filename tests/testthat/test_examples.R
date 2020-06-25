context("Generic unit testing")
library(MiRNAQCD)

datasetAlpha <- read.table(file="example_dataset_alpha.dat", header=TRUE)
preprocessedDatasetAlpha <- miRNA_expressionPreprocessing(datasetAlpha, multipletSize=3)
filename_qt=tempfile()
filename_threshold=tempfile()

test_that("Test quality thresholds", {
	expect_equal(length(preprocessedDatasetAlpha[,1]), 357)
	qualityThresholdAlpha <- miRNA_assessQualityThreshold(preprocessedDatasetAlpha, significanceLevel=0.05, saveOutputFile=TRUE, outputFileName=filename_qt)
	expect_equal(as.character(qualityThresholdAlpha[[2]][1]), "0.51")
	expect_equal(as.character(qualityThresholdAlpha[[2]][2]), "0.5")
	expect_equal(as.character(qualityThresholdAlpha[[2]][3]), "0.53")
})

test_that("Test load QT & outlier removal", {
	qualityThresholdAlpha <- miRNA_loadQualityThreshold(filename_qt)
 	purifiedDatasetAlpha <- miRNA_removeOutliers(preprocessedDatasetAlpha, qualityThresholdAlpha)
	expect_equal(length(purifiedDatasetAlpha[,1]), 345)
})

qualityThresholdAlpha <- miRNA_loadQualityThreshold(filename_qt)
purifiedDatasetAlpha <- miRNA_removeOutliers(preprocessedDatasetAlpha, qualityThresholdAlpha)
Target <- c("A")
Versus <- c("B", "C")

test_that("Test analysis mode", {
	statisticsAlpha <- miRNA_classifierSetup(purifiedDatasetAlpha, inputTargetList=Target, inputVersusList=Versus, saveOutputFile=FALSE)
	expect_equal(as.character(statisticsAlpha[1,1]), "FX")
	expect_equal(as.character(statisticsAlpha[1,3]), "55")
	expect_equal(as.character(statisticsAlpha[1,6]), "0.17")
})

mirnaToUse <- c("FX", "FZ")
coefficientsToUse <- c(1.0, -1.0)

test_that("Test training mode", {
	thresholdValues <- miRNA_classifierSetup(purifiedDatasetAlpha, inputTargetList=Target, inputVersusList=Versus, inputMiRNAList=mirnaToUse, coeffList=coefficientsToUse, saveOutputFile=TRUE, outputFileBasename=filename_threshold)
	expect_equal(thresholdValues[1,1], 4.55, tolerance=0.02)
	expect_equal(thresholdValues[1,4], 0.95, tolerance=0.02)
})

datasetBeta <- read.table(file="example_dataset_beta.dat", header=TRUE)
preprocessedDatasetBeta <- miRNA_expressionPreprocessing(datasetBeta, multipletSize=3)
purifiedDatasetBeta <- miRNA_removeOutliers(preprocessedDatasetBeta, qualityThresholdAlpha)
thresholdValues <- miRNA_loadDiagnosticThreshold(paste(filename_threshold, ".txt", sep=""))

test_that("Test check dataset beta", {
	expect_equal(length(preprocessedDatasetBeta[,1]), 597)
	expect_equal(length(purifiedDatasetBeta[,1]), 568)
})

test_that("Test diagnosis", {
	diagnosedBeta <- miRNA_diagnosis(purifiedDatasetBeta, inputMiRNAList=mirnaToUse, coeffList=coefficientsToUse, inputThreshold=thresholdValues, saveOutputFile=FALSE)
	expect_equal(diagnosedBeta[1,3], 0.456)
	expect_equal(diagnosedBeta[100,3], 4.7996)
	system("rm -f .temp*")
})
