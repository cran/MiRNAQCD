## --------------------------------------------------------------------------
##
## This file is part of the miRNA-QC-and-Diagnosis software package.
##
## Version 1.0 - June 2020
##
##
## The miRNA-QC-and-Diagnosis package is free software; you can use it,
## redistribute it, and/or modify it under the terms of the GNU General
## Public License version 3 as published by the Free Software Foundation.
## The full text of the license can be found in the file LICENSE.txt at the top
## level of the package distribution.
##
## Authors:
##	Michele Castelluzzo (1), Alessio Perinelli (1), Simone Detassis (2),
##	Michela A. Denti (2) and Leonardo Ricci (1,3)
##	(1) Department of Physics, University of Trento, 38123 Trento, Italy
##	(2) Department of Cellular, Computational and Integrative Biology
##		(CIBIO), University of Trento, 38123 Trento, Italy
##	(3) CIMeC, Center for Mind/Brain Sciences, University of Trento,
##		38068 Rovereto, Italy
##
##	michele.castelluzzo@unitn.it
##	alessio.perinelli@unitn.it
##	michela.denti@unitn.it
##	leonardo.ricci@unitn.it
##	https://github.com/LeonardoRicci/
##	https://nse.physics.unitn.it/
##
##
## If you use the miRNA-QC-and-Diagnosis package for your analyses, please cite:
##
##	L. Ricci, V. Del Vescovo, C. Cantaloni, M. Grasso, M. Barbareschi and
##	M. A. Denti, Statistical analysis of a Bayesian classifier based on the
##	expression of miRNAs, BMC Bioinformatics 16:287 (2015).
##	DOI: 10.1186/s12859-015-0715-9
##
##
## --------------------------------------------------------------------------

#' Removal of dataset outliers.
#'
#' This function removes outliers from a given dataset according to a set of quality threshold values.
#'
#' @param inputDataset Dataset (data frame) to be cleaned of outliers. The data frame must comply with the output format of the preprocessing function (miRNA_expressionPreprocessing), thus containing the columns 'Subject', 'miRNA', 'Mean', 'StdDev', 'SampleSize' and possibly 'Class'. Any other column is ignored, and any missing column forbids execution.
#' @param qualityThresholdFrame Critical sigma values (data frame) to be used. The data frame must comply with the output format of the ebbc function for critical sigma assessment (miRNA_assessQualityThreshold), thus containing the columns 'miRNA' and 'QualityThreshold'. Any other column is ignored, and any missing column forbids execution.
#'
#' Beware! Entries of the dataset for which 'miRNA' is not present in the data frame of critical sigma values are copied in output without any filtering.
#'
#' @return A data frame corresponding to a copy of the input dataset devoid of outliers. The output data frame thus contains the columns 'Subject', 'miRNA', 'Mean', 'StdDev', 'Variance', 'SampleSize' and possibly 'Class'.
#'
#' @examples
#' requiredDataFile = paste(system.file(package="MiRNAQCD"),
#'			"/extdata/dataset_alpha_prep.dat", sep='')
#' myDataFrame <- read.table(file=requiredDataFile, header=TRUE)
#' requiredQtFile = paste(system.file(package="MiRNAQCD"),
#'			"/extdata/dataset_alpha_qt.dat", sep='')
#' qtDataFrame <- read.table(file=requiredQtFile, header=TRUE)
#' myDataFrameCleaned <- miRNA_removeOutliers(myDataFrame, qtDataFrame)

#' @export
miRNA_removeOutliers <- function(inputDataset, qualityThresholdFrame) {

	if (!(("Subject" %in% colnames(inputDataset)) & ("miRNA" %in% colnames(inputDataset)) & ("Mean" %in% colnames(inputDataset)) & ("StdDev" %in% colnames(inputDataset)) & ("SampleSize" %in% colnames(inputDataset))))  {
		stop("ERROR: unsuitable dataset format. Dataset must contain columns 'Subject', 'miRNA', 'Mean', 'StdDev', 'SampleSize'.\n")
	}

	if (!(("miRNA" %in% colnames(qualityThresholdFrame)) & ("QualityThreshold" %in% colnames(qualityThresholdFrame))))  {
		stop("ERROR: unsuitable qualityThresholdFrame format. Data frame must contain columns 'miRNA', 'QualityThreshold'.\n")
	}

	if (length(inputDataset[1,]) > 7) {
		warning("WARNING: more than 7 dataset columns. Columns other than 'Subject', 'miRNA', 'Mean', 'StdDev', 'SampleSize', 'Class' will not be present in output.\n")
	}

	if ("Class" %in% colnames(inputDataset)) {
		selectedDataset <- subset(inputDataset, select=c("Subject", "miRNA", "Mean", "StdDev", "SampleSize", "Class"))
	} else {
		selectedDataset <- subset(inputDataset, select=c("Subject", "miRNA", "Mean", "StdDev", "SampleSize"))
	}

	setOfFeatures <- unique(qualityThresholdFrame$miRNA)
	setOfFeatures <- setOfFeatures[order(setOfFeatures)]
	selectedDataset <- selectedDataset[selectedDataset$miRNA %in% setOfFeatures, ]

	for(i in seq (1,nrow(qualityThresholdFrame))) {
		feature <- qualityThresholdFrame$miRNA[i]
		sigmaLimit <- as.numeric(as.vector(qualityThresholdFrame$QualityThreshold[i]))
		selectedDataset <- selectedDataset[!(selectedDataset$miRNA == feature & selectedDataset$StdDev > sigmaLimit),]
	}

	return(selectedDataset)
}
